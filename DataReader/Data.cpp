#include "DataReader/Data.hpp"
Data::Data(bool binning, unsigned int maxBins,double maxW) :
fBinned(binning), fmaxBins(maxBins), maxWeight(maxW)
{

}

void rndReduceSet(unsigned int size, std::shared_ptr<Generator> gen,
		Data* in1, Data* out1, Data*  in2, Data*  out2){
	if(!in1 )
		throw std::runtime_error("rndSubSet() | No input data set!");
	if(!out1 )
		throw std::runtime_error("rndSubSet() | No output data set!");
	if( out1->getNEvents() )
		throw std::runtime_error("rndSubSet() | First output sample not empty!");
	if(in2){
		if(in1->getNEvents() != in2->getNEvents())
			throw std::runtime_error("rndSubSet() | Samples have different event count!");
		if(!out2)
			throw std::runtime_error("rndSubSet() | Second output set is NULL!");
		if( out2->getNEvents() )
			throw std::runtime_error("rndSubSet() | Second output sample not empty!");
	}

	unsigned int totalSize = in1->getNEvents();
	unsigned int newSize = totalSize;
	/* 1th method: new Sample has exact size, but possibly events are added twice.
	 * We would have to store all used events in a vector and search the vector at every event -> slow */
	/*unsigned int t=0;
	 unsigned int d=0;
	 while(t<newSize){
	 d = (unsigned int) gen->getUniform()*totalSize;
	 newSample->pushEvent(fEvents[d]);
	 t++;
	 }*/

	/* 2nd method: events are added once only, but total size of sample varies with sqrt(N) */
	for(unsigned int i=0; i<totalSize; i++){//count how many events are not within PHSP
		dataPoint point;
		try{
			point = dataPoint(in1->getEvent(i));
		} catch (BeyondPhsp& ex){ //event outside phase, remove
			newSize--;
			continue;
		}
		//		dataPoint point(in1->getEvent(i));
		//		if(!Kinematics::instance()->isWithinPhsp(point)) newSize--;
	}
	double threshold = (double)size/newSize; //calculate threshold
	for(unsigned int i=0; i<totalSize; i++){
		dataPoint point;
		try{
			point = dataPoint(in1->getEvent(i));
		} catch (BeyondPhsp& ex){ //event outside phase, remove
			continue;
		}
		//		dataPoint point(in1->getEvent(i)); //use first sample for hit&miss
		//		if(!Kinematics::instance()->isWithinPhsp(point)) continue;
		if( gen->getUniform() < threshold ){
			out1->pushEvent(in1->getEvent(i));
			//write second sample if event from first sample were accepted
			if(in2) out2->pushEvent(in2->getEvent(i));
		}
	}
	if(out2) assert(out1->getNEvents() == out2->getNEvents());

	BOOST_LOG_TRIVIAL(debug) << "DataReader::rndReduceSet() | sample size reduced to "<<out1->getNEvents();
	return;
}

std::shared_ptr<Data> Data::rndSubSet(unsigned int size, std::shared_ptr<Generator> gen){
	std::shared_ptr<Data> out(this->EmptyClone());
	rndReduceSet(size, gen,this, out.get());
	return out;
}



void Data::resetWeights(double w)
{
	for(unsigned int i=0; i<fEvents.size(); i++)
		fEvents.at(i).setWeight(w);
	maxWeight=w;
	return;
}

double Data::getMaxWeight() const
{
	return maxWeight;
}

void Data::reduceToPhsp()
{
	std::vector<Event> tmp;
	BOOST_LOG_TRIVIAL(info)<<"Data::reduceToPhsp() | "
			"Remove all events outside PHSP boundary from data sample.";

	for(unsigned int evt=0; evt<fEvents.size(); evt++){
		dataPoint point;
		try{
			point = dataPoint(fEvents.at(evt));
		} catch (BeyondPhsp& ex){ //event outside phase, remove
			continue;
		}
		//		if(Kinematics::instance()->isWithinPhsp(p))
		tmp.push_back(fEvents.at(evt));
	}
	BOOST_LOG_TRIVIAL(info)<<"Data::reduceToPhsp() | "
			<<tmp.size()<<" from "<<fEvents.size()
			<<"("<<((double)tmp.size())/fEvents.size()*100 <<"%) were kept.";
	fEvents = tmp;
	return;
}
void Data::resetEfficiency(double e)
{
	for(unsigned int evt=0; evt<fEvents.size(); evt++){
		fEvents.at(evt).setEfficiency(e);
	}
}
void Data::reduce(unsigned int newSize)
{
	if(newSize >= fEvents.size()) {
		BOOST_LOG_TRIVIAL(error) << "RooReader::reduce() requested size too large, cant reduce sample!";
		return;
	}
	fEvents.resize(newSize);
}

void Data::setEfficiency(std::shared_ptr<Efficiency> eff)
{
	for(unsigned int evt=0; evt<fEvents.size(); evt++){
		dataPoint point;
		try{
			point = dataPoint(fEvents.at(evt));
		} catch (BeyondPhsp& ex){ //event outside phase, remove
			continue;
		}
		//		dataPoint point(fEvents.at(evt));
		double val = eff->evaluate(point);
		fEvents.at(evt).setEfficiency(val);
	}
}
void Data::Clear()
{
	fEvents.clear();
}

bool Data::hasWeights()
{
	bool has=0;
	for(unsigned int evt=0; evt<fEvents.size(); evt++){
		if(fEvents.at(evt).getWeight()!=1.) {
			has=1;
			break;
		}
	}
	return has;
}

const ParameterList& Data::getListOfData()
{
	//dataList already filled? return filled one
	if(dataList.GetNParameter() != 0)
		return dataList;

	int size = Kinematics::instance()->GetNVars();
	std::vector< std::vector<double> > data(size, std::vector<double>());
	std::vector< double > eff;
	eff.reserve(fEvents.size());
	std::vector< double > weight;
	weight.reserve(fEvents.size());

	auto itr = fEvents.begin();
	for( ; itr!=fEvents.end(); ++itr){
		dataPoint point(*itr);
		eff.push_back(point.getEfficiency());
		weight.push_back(point.getWeight());
		for(int i=0; i<size; ++i)
			data.at(i).push_back(point.getVal(i));
	}

	//Add data vector to ParameterList
	for(int i=0; i<size; ++i){
		std::shared_ptr<MultiDouble> tmp(
				new MultiDouble(
						Kinematics::instance()->getVarNames().at(i),
						data.at(i)
				)
		);
		dataList.AddParameter(tmp);
	}

	//Adding efficiency at the end
	dataList.AddParameter(
			std::shared_ptr<MultiDouble>(
					new MultiDouble(
							"Efficiency",
							eff
					)
			)
	);

	//Adding weight at the end
	dataList.AddParameter(
			std::shared_ptr<MultiDouble>(
					new MultiDouble(
							"Weight",
							weight
					)
			)
	);

	return dataList;
}

std::vector<dataPoint> Data::getDataPoints() const
{
	std::vector<dataPoint> vecPoint;
	for(int i=0; i<fEvents.size(); i++){
		dataPoint point;
		try{
			point = dataPoint(fEvents.at(i));
		} catch (BeyondPhsp& ex){ //event outside phase, remove
			continue;
		}
		vecPoint.push_back(point);
	}
	return vecPoint;
}

void Data::setResolution(std::shared_ptr<Resolution> res)
{
	for(int i=0; i<fEvents.size(); i++)
		res->resolution(fEvents.at(i));
}

void Data::Add(Data& otherSample)
{
	std::vector<Event> otherEvents = otherSample.getEvents();
	fEvents.insert(fEvents.end(), otherEvents.begin(), otherEvents.end());
	if(otherSample.getMaxWeight() > maxWeight)
		maxWeight = otherSample.getMaxWeight();
	return;
}

void Data::applyCorrection(DataCorrection& corr)
{
	double sumWeightSq=0;
	for(int i=0; i<fEvents.size(); i++){
		double w = corr.getCorrection(fEvents.at(i));
		if( w < 0 )
			throw std::runtime_error("Data::applyCorrection() | "
					"Negative weight!");
		sumWeightSq += w*w;
		double oldW = fEvents.at(i).getWeight();
		if(w*oldW > maxWeight) maxWeight = w*oldW;
		fEvents.at(i).setWeight(w*oldW);
	}
	BOOST_LOG_TRIVIAL(info)<<"Data::applyCorrection() | "
			"Sample corrected! Sum of weights squared is "<<sumWeightSq;
	return;
}

const int Data::getBin(const int i, double& m12, double& weight)
{
	if(!fBinned) return -1;

	m12 = fBins[i].first;
	weight = fBins[i].second;

	return 1;
}

//! Add event to data sample
void Data::pushEvent(const Event& evt)
{
	fEvents.push_back(evt);
	if( evt.getWeight() > maxWeight ) maxWeight = evt.getWeight();
}
