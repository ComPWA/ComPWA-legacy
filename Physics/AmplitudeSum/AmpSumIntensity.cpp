//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding flatte type resonance, removing root dependence
//-------------------------------------------------------------------------------

#include "Core/PhysConst.hpp"
#include "Core/Functions.hpp"

// Boost header files go here
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Physics/AmplitudeSum/NonResonant.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

#include "gsl/gsl_monte_vegas.h"

#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"

AmpSumIntensity::AmpSumIntensity(normStyle ns, std::shared_ptr<Efficiency> eff,
		unsigned int nCalls) :
		_normStyle(ns), _calcMaxFcnVal(0), eff_(eff), _nCalls(nCalls)
{
	result.AddParameter(
			std::shared_ptr<DoubleParameter>(
					new DoubleParameter("AmpSumResult")
			)
	);
	return;
}

//! Configure resonance from ptree
void AmpSumIntensity::Configure(const boost::property_tree::ptree &pt)
{
	BOOST_FOREACH( ptree::value_type const& v,
			pt.get_child("amplitude_setup") ) {
		/* We try to configure each type of resonance. In case that v.first does
		 * not contain the correct string, a BadConfig is thrown. We catch it
		 * and try the next type.
		 */
		try{
			auto amp=new AmpRelBreitWignerRes(_normStyle, _nCalls);
			amp->Configure(v,params);
			resoList.push_back(std::shared_ptr<AmpAbsDynamicalFunction>(amp));
			BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::Configure() | "
					"adding amplitude:"
					<<std::endl<<amp->to_str();
		} catch (BadConfig& ex) {}
		try{
			auto amp=new AmpFlatteRes(_normStyle, _nCalls);
			amp->Configure(v,params);
			resoList.push_back(std::shared_ptr<AmpAbsDynamicalFunction>(amp));
			BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::Configure() | "
					"adding amplitude:"
					<<std::endl<<amp->to_str();
		} catch (BadConfig& ex) {}
		try{
			auto amp=new NonResonant(_normStyle, _nCalls);
			amp->Configure(v,params);
			resoList.push_back(std::shared_ptr<AmpAbsDynamicalFunction>(amp));
			BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::Configure() | "
					"adding amplitude:"
					<<std::endl<<amp->to_str();
		} catch (BadConfig& ex) {}
	}
	BOOST_LOG_TRIVIAL(info) << "AmpSumIntensity::Configure() | "
			"Loaded property tree!";
	return;
}

//! Save resonance from to ptree
void AmpSumIntensity::Save(std::string filePath)
{
	using boost::property_tree::ptree;
	boost::property_tree::ptree ptSub, pt;

	BOOST_LOG_TRIVIAL(info) << "AmpSumIntensity::Save() | "
			"Save amplitude to "<<filePath<<"!";
	auto it = resoList.begin();
	for( ; it!= resoList.end(); ++it)
		(dynamic_pointer_cast<AmpAbsDynamicalFunction>(*it))->Save(ptSub);

	pt.add_child("amplitude_setup",ptSub);
	//new line at the end
	boost::property_tree::xml_writer_settings<char> settings('\t', 1);
	// Write the property tree to the XML file.
	write_xml(filePath, pt,std::locale(), settings);

}

//! Getter function for function tree
virtual std::shared_ptr<FunctionTree> AmpSumIntensity::GetTree(
		ParameterList& sample, ParameterList& toySample,
		std::string suffix="")
{
	return setupBasicTree(sample,toySample, suffix);
}
		//we assume that the total efficiency of the sample is stored as efficiency of each event
//		accSampleEff = mAccSample.GetMultiDouble(effId)->GetValue(0);
//		BOOST_LOG_TRIVIAL(info)<<"MinLogLH::MinLogLH() total efficiency of unbinned correction sample is set to "<<accSampleEff;

	//Store position of efficiency and weight vector
//	effId = Kinematics::instance()->GetNVars();
//	weightId = Kinematics::instance()->GetNVars()+1;

//void MinLogLH::iniLHtree(){
//	//reset all trees before generating new trees; saves a lot of virtual memory
//	signalPhspTree = std::shared_ptr<FunctionTree>();
//	signalTree_amp = std::shared_ptr<FunctionTree>();
//	signalPhspTree_amp = std::shared_ptr<FunctionTree>();
//	bkgPhspTree = std::shared_ptr<FunctionTree>();
//	bkgTree_amp = std::shared_ptr<FunctionTree>();
//	bkgPhspTree_amp = std::shared_ptr<FunctionTree>();
//
//	BOOST_LOG_TRIVIAL(debug) << "MinLogLH::iniLHtree() constructing the LH tree";
//
//	if(useFunctionTree) return;
//	if(!amp->hasTree()){
//		throw std::runtime_error("MinLogLH::iniLHtree() amplitude has no tree");
//	}
//	if(ampBkg && !ampBkg->hasTree()){
//		throw std::runtime_error("MinLogLH::iniLHtree() amplitude has no tree");
//	}
//
//	//----Strategies needed
//	std::shared_ptr<Strategy> mmultStrat(new MultAll(ParType::MCOMPLEX));
//	std::shared_ptr<Strategy> mmultDStrat(new MultAll(ParType::MDOUBLE));
//	std::shared_ptr<Strategy> multiDoubleAddStrat(new AddAll(ParType::MDOUBLE));
//	std::shared_ptr<Strategy> multiComplexAddStrat(new AddAll(ParType::MCOMPLEX));
//	std::shared_ptr<Strategy> msqStrat(new AbsSquare(ParType::MDOUBLE));
//	std::shared_ptr<Strategy> mlogStrat(new LogOf(ParType::MDOUBLE));
//	std::shared_ptr<Strategy> multStrat(new MultAll(ParType::COMPLEX));
//	std::shared_ptr<Strategy> multDStrat(new MultAll(ParType::DOUBLE));
//	std::shared_ptr<Strategy> addStrat(new AddAll(ParType::DOUBLE));
//	std::shared_ptr<Strategy> addComplexStrat(new AddAll(ParType::COMPLEX));
//	std::shared_ptr<Strategy> sqStrat(new AbsSquare(ParType::DOUBLE));
//	std::shared_ptr<Strategy> logStrat(new LogOf(ParType::DOUBLE));
//	std::shared_ptr<Strategy> complStrat(new Complexify(ParType::COMPLEX));
//	std::shared_ptr<Strategy> invStrat(new Inverse(ParType::DOUBLE));
//
//	BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() construction normalization tree";
//
//	//=== Signal normalization
//	signalPhspTree = std::shared_ptr<FunctionTree>(new FunctionTree());
//	signalPhspTree->createHead("invNormLH", invStrat);// 1/normLH
//	// normLH = phspVolume/N_{mc} |T_{evPHSP}|^2
//	signalPhspTree->createNode("normFactor", multDStrat, "invNormLH");
//	// sumAmp = \sum_{evPHSP} |T_{evPHSP}|^2
//	signalPhspTree->createNode("sumAmp", addStrat,"normFactor");
//	std::shared_ptr<MultiDouble> eff, weightPhsp;
//	signalPhspTree->createLeaf("phspVolume",
//			Kinematics::instance()->getPhspVolume(), "normFactor");
//
//	//Which kind of efficiency correction should be used?
//	if(!accSample) {//binned
//		int phspSampleSize = mPhspSample.GetMultiDouble(0)->GetNValues();
//		signalPhspTree_amp = amp->GetTree(mPhspSample,mPhspSample,"_Phsp");
//
//		weightPhsp = mPhspSample.GetMultiDouble(weightId);
//		double sumWeights =
//				std::accumulate(weightPhsp->Begin(), weightPhsp->End(), 0.0);
//		signalPhspTree->createLeaf("InvNmc",
//				1/ ( (double) sumWeights ), "normFactor");
//		signalPhspTree->createNode("IntensPhspEff", mmultDStrat, "sumAmp",
//				phspSampleSize, false); //|T_{ev}|^2
//		eff = mPhspSample.GetMultiDouble(effId);
//		signalPhspTree->createLeaf("eff", eff, "IntensPhspEff"); //efficiency
//		signalPhspTree->createLeaf("weightPhsp", weightPhsp, "IntensPhspEff");
//		signalPhspTree->createNode("IntensPhsp", msqStrat, "IntensPhspEff",
//				phspSampleSize, false); //|T_{ev}|^2
//		BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() | Setting up "
//				"normalization tree, using toy sample and assume that "
//				"efficiency values are saved for every event!";
//
//		signalPhspTree->insertTree(signalPhspTree_amp, "IntensPhsp");
//	}
//	else {//unbinned
//		int accSampleSize = mAccSample.GetMultiDouble(0)->GetNValues();
//		signalPhspTree->createNode("weightIntensPhsp", mmultDStrat, "sumAmp",
//				accSampleSize, false);
//		weightPhsp = mPhspSample.GetMultiDouble(weightId);
//		signalPhspTree->createLeaf("weightPhsp", weightPhsp, "weightIntensPhsp");
//		signalPhspTree->createNode("IntensPhsp", msqStrat, "weightIntensPhsp",
//				accSampleSize, false); //|T_{ev}|^2
//		double sumWeights =
//				std::accumulate(weightPhsp->Begin(), weightPhsp->End(), 0.0);
//		signalPhspTree->createLeaf("InvNmc",
//				1/ ( (double) sumWeights/accSampleEff ), "normFactor");
//		signalPhspTree_amp = amp->GetTree(mAccSample,mPhspSample,"_Phsp");
//		BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() setting up normalization tree, "
//				"using sample of accepted phsp events for efficiency correction!";
//		signalPhspTree->insertTree(signalPhspTree_amp, "IntensPhsp"); //Sum of resonances, at each point
//	}
//	//=== Background normalization
//	if(ampBkg){
//		bkgPhspTree = std::shared_ptr<FunctionTree>(new FunctionTree());
//		bkgPhspTree->createHead("invBkgNormLH", invStrat);// 1/normLH
//		bkgPhspTree->createNode("normFactor", multDStrat, "invBkgNormLH"); // normLH = phspVolume/N_{mc} |T_{evPHSP}|^2
//		bkgPhspTree->createNode("sumAmp", addStrat,"normFactor"); // sumAmp = \sum_{evPHSP} |T_{evPHSP}|^2
//		bkgPhspTree->createLeaf("phspVolume", Kinematics::instance()->getPhspVolume(), "normFactor");
//		if(!accSample) {//binned
//			int phspSampleSize = mPhspSample.GetMultiDouble(0)->GetNValues();
//			bkgPhspTree_amp = ampBkg->GetTree(mPhspSample,mPhspSample,"_Phsp");
//
//			weightPhsp = mPhspSample.GetMultiDouble(weightId);
//			double sumWeights =
//					std::accumulate(weightPhsp->Begin(), weightPhsp->End(), 0.0);
//			bkgPhspTree->createLeaf("InvNmc", 1/ ( (double) sumWeights), "normFactor");
//			bkgPhspTree->createNode("IntensPhspEff", mmultDStrat, "sumAmp",
//					phspSampleSize, false); //|T_{ev}|^2
//			eff = mPhspSample.GetMultiDouble(effId);
//			bkgPhspTree->createLeaf("eff", eff, "IntensPhspEff"); //efficiency
//			bkgPhspTree->createLeaf("weightPhsp", weightPhsp, "IntensPhspEff");
//			bkgPhspTree->createNode("IntensPhsp", msqStrat, "IntensPhspEff",
//					phspSampleSize, false); //|T_{ev}|^2
//			BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() | Setting up tree "
//					"for background normalization, using toy sample and assume "
//					"that efficiency values are saved for every event!";
//
//			//Sum of resonances, at each point
//			bkgPhspTree->insertTree(bkgPhspTree_amp, "IntensPhsp");
//		}
//		else {//unbinned
//			int accSampleSize = mAccSample.GetMultiDouble(0)->GetNValues();
//			bkgPhspTree->createNode("weightIntensPhsp", mmultDStrat, "sumAmp",
//					accSampleSize, false);
//			weightPhsp = mPhspSample.GetMultiDouble(weightId);
//			bkgPhspTree->createLeaf("weightPhsp", weightPhsp, "weightIntensPhsp");
//			bkgPhspTree->createNode("IntensPhsp", msqStrat, "weightIntensPhsp",
//					accSampleSize, false); //|T_{ev}|^2
//			double sumWeights =
//					std::accumulate(weightPhsp->Begin(), weightPhsp->End(), 0.0);
//			bkgPhspTree->createLeaf("InvNmc",
//					1/ ( (double)sumWeights/accSampleEff ), "normFactor");
//			bkgPhspTree_amp = ampBkg->GetTree(mAccSample,mPhspSample,"_Phsp");
//			BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() | Setting up tree "
//					"for background normalization, using sample of accepted "
//					"phsp events for efficiency correction!";
//
//			bkgPhspTree->insertTree(bkgPhspTree_amp, "IntensPhsp"); //Sum of resonances, at each point
//		}
//	}
//
//	BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() construction LH tree";
//	/* CONSTRUCTION OF THE LIKELIHOOD:
//	 * We denote the coherent sum over all resonances with T:
//	 * 		T := \sum_{i,j} c_i c_j^*A_iA_j^*
//	 * The negative log LH is given by:
//	 * 		-log L = - N/(\sum_{ev} w_{ev}) \sum_{ev} w_{ev} \log{f_{bkg} \frac{|T|^2}{\int_{DP} |T|^2} + (1-f_{bkg})}
//	 * The sum over all weights is necessary to normalize the weights to one. Otherwise the error
//	 * estimate is incorrect. The LH normalization is norm_{LH} = \int_{DP} |T|^2.
//	 * This formulation includes event weights as well as a flat background description. f_{bkg} is
//	 * the fraction of background in the sample. Using both is of course non-sense. Set weights to
//	 * one OR f_{bkg} to zero.
//	 */
//	physicsTree = std::shared_ptr<FunctionTree>(new FunctionTree());
//	int sampleSize = mData.GetMultiDouble(0)->GetNValues();
//	/* Setup basic tree
//	 * head node is 'Amplitude' which contains the complex amplitude values for each event in sample
//	 */
//	signalTree_amp = amp->GetTree(mData,mPhspSample,"data");
//	std::shared_ptr<MultiDouble> weight = mData.GetMultiDouble(weightId);
//
//	//-log L = (-1)*N/(\sum_{ev} w_{ev}) \sum_{ev} ...
//	physicsTree->createHead("LH", multDStrat);
//	physicsTree->createLeaf("minusOne", -1 ,"LH");
//	physicsTree->createLeaf("nEvents", sampleSize ,"LH");
//	physicsTree->createNode("invSumWeights", invStrat,"LH");
//	physicsTree->createNode("sumEvents", addStrat, "LH");
//	physicsTree->createNode("sumWeights", addStrat, "invSumWeights");
//	physicsTree->createLeaf("weight", weight, "sumWeights");
//	physicsTree->createNode("weightLog", mmultDStrat, "sumEvents",
//			sampleSize, false); //w_{ev} * log( I_{ev} )
//	physicsTree->createLeaf("weight", weight, "weightLog");
//	physicsTree->createNode("Log", mlogStrat, "weightLog", sampleSize, false);
//	//I_{ev} = x_{ev} + (1-f_{bkg})
//	physicsTree->createNode("Add", multiDoubleAddStrat, "Log",
//			sampleSize, false);
//
//	//signal term
//	physicsTree->createNode("normIntens", mmultDStrat, "Add", sampleSize, false);// x=f_{bkg}|T|^2/norm_{LH}
//	physicsTree->createLeaf("signalFrac", signalFraction, "normIntens");
//	physicsTree->insertTree(signalPhspTree, "normIntens"); //provides 1/normLH
//	physicsTree->createNode("Intens", msqStrat, "normIntens", sampleSize, false);
//	physicsTree->insertTree(signalTree_amp,"Intens");
//	//background term
//	if(ampBkg){
//		physicsTree->createNode("normBkg", mmultDStrat, "Add", sampleSize, false);// x=f_{bkg}|T|^2/norm_{LH}
//		physicsTree->createLeaf("OneMinusBkgFrac", (1-signalFraction), "normBkg");
//		bkgTree_amp= ampBkg->GetTree(mData,mPhspSample,"data");
//		physicsTree->insertTree(bkgPhspTree, "normBkg"); //provides 1/normLH
//		physicsTree->createNode("IntensBkg", msqStrat, "normBkg", sampleSize, false);
//		physicsTree->insertTree(bkgTree_amp,"IntensBkg");
//	}
//	physicsTree->recalculate();
//	//	std::string treeString = physicsTree->head()->to_str(10);
//	//	BOOST_LOG_TRIVIAL(debug) << std::endl << treeString;
//	if(!physicsTree->sanityCheck()) {
//		throw std::runtime_error("MinLogLH::iniLHtree() tree has structural problems. Sanity check not passed!");
//	}
//	BOOST_LOG_TRIVIAL(debug) <<"MinLogLH::iniLHtree() construction of LH tree finished!";
//	useFunctionTree=1;
//	return;
//}
std::shared_ptr<FunctionTree> AmpSumIntensity::setupBasicTree(
		ParameterList& sample, ParameterList& toySample,std::string suffix)
{
	int sampleSize = sample.GetMultiDouble(0)->GetNValues();
	int toySampleSize = toySample.GetMultiDouble(0)->GetNValues();

	BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::setupBasicTree() generating new tree!";
	if(sampleSize==0){
		BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::setupBasicTree() data sample empty!";
		return std::shared_ptr<FunctionTree>();
	}
	if(toySampleSize==0){
		BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::setupBasicTree() toy sample empty!";
		return std::shared_ptr<FunctionTree>();
	}

	//------------Setup Tree---------------------
	std::shared_ptr<FunctionTree> newTree(new FunctionTree());

	//----Strategies needed
	std::shared_ptr<AddAll> maddStrat(new AddAll(ParType::MCOMPLEX));
	newTree->createHead("Amplitude"+suffix, maddStrat, sampleSize);

	auto it = resoList.begin();
	for( ; it!=resoList.end(); ++it){
		if(!(*it)->GetEnable()) continue;
		std::shared_ptr<FunctionTree> resTree= (*it)->SetupTree(sample,
				toySample, ""+(*it)->GetName());
		if(!resTree->sanityCheck())
			throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
					"Resonance tree didn't pass sanity check!");
		resTree->recalculate();
		newTree->insertTree(resTree, "Amplitude"+suffix);
	}

	BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::setupBasicTree(): tree constructed!!";
	return newTree;
}

double AmpSumIntensity::getMaxVal(std::shared_ptr<Generator> gen)
{
	if(!_calcMaxFcnVal) calcMaxVal(gen);
	return _maxFcnVal;
}

double AmpSumIntensity::getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen)
{
	calcMaxVal(par,gen);
	return _maxFcnVal;
}

void AmpSumIntensity::calcMaxVal(ParameterList& par, std::shared_ptr<Generator> gen)
{
	setParameterList(par);
	return calcMaxVal(gen);
}

void AmpSumIntensity::calcMaxVal(std::shared_ptr<Generator> gen)
{
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());

	double maxM23=-999; double maxM13=-999; double maxVal=0;
	for(unsigned int i=0; i<_nCalls; i++){
		auto m13sq_limit = kin->GetMinMax(1);
		auto m23sq_limit = kin->GetMinMax(0);

		double m23sq = gen->getUniform()*(m23sq_limit.second-m23sq_limit.first)
										+m23sq_limit.first;
		double m13sq = gen->getUniform()*(m13sq_limit.second-m13sq_limit.first)
										+m13sq_limit.first;
		dataPoint point;
		try{
			Kinematics::instance()->FillDataPoint(1,0,m13sq,m23sq,point);
		} catch (BeyondPhsp& ex){
			if(i>0) i--;
			continue;
		}
		ParameterList res = intensity(point);
		double intens = *res.GetDoubleParameter(0);
		if(intens>maxVal) {
			maxM23=m23sq; maxM13=m13sq;
			maxVal=intens;
		}
	}
	_maxFcnVal=maxVal;
	_calcMaxFcnVal=1;
	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity::calcMaxVal() calculated maximum of amplitude: "
			<<_maxFcnVal<<" at m23sq="<<maxM23<<"/m13sq="<<maxM13;
	return ;
}

double AmpSumIntensity::evaluate(double x[], size_t dim)
{
	/* Calculation amplitude integral (excluding efficiency) */
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//set data point: we assume that x[0]=m13 and x[1]=m23
	dataPoint point;
	try{
		kin->FillDataPoint(1,0,x[0],x[1],point);
	} catch (BeyondPhsp& ex){
		return 0;
	}
	ParameterList res = intensityNoEff(point);
	double intens = *res.GetDoubleParameter(0);
	return intens;
}

double evalWrapperAmpSumIntensity(double* x, size_t dim, void* param)
{
	/* We need a wrapper here because intensity() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	return static_cast<AmpSumIntensity*>(param)->evaluate(x,dim);
}

const double AmpSumIntensity::integral(ParameterList& par)
{
	setParameterList(par);
	return integral();
}

const double AmpSumIntensity::integral()
{
	/* Integration functionality was tested with a model with only one normalized amplitude.
	 * The integration result is equal to the amplitude coefficient^2.
	 */
	size_t dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	auto m13sq_limit = kin->GetMinMax(1);
	auto m23sq_limit = kin->GetMinMax(0);
	double xLimit_low[2] = {m13sq_limit.first,m23sq_limit.first};
	double xLimit_high[2] = {m13sq_limit.second,m23sq_limit.second};
	//	double xLimit_low[2] = {0,0};
	//	double xLimit_high[2] = {10,10};
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function G = {&evalWrapperAmpSumIntensity,dim, const_cast<AmpSumIntensity*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&G, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::integrate() Integration result for amplitude sum: "
			<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

double interferenceIntegralWrapper(double* x, size_t dim, void* param)
{
	/* Calculation amplitude integral (excluding efficiency) */
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	dataPoint point;
	try{
		Kinematics::instance()->FillDataPoint(0,1,x[1],x[0],point);
	} catch (BeyondPhsp& ex){
		return 0;
	}

	AmpSumIntensity* amp = static_cast<AmpSumIntensity*>(param);

	auto A = amp->tmpA;
	auto B = amp->tmpB;
	double intens = (
			(*A)->Evaluate(point)*std::conj((*B)->Evaluate(point))
	).real();
	if( A == B ) return intens;
	return 2*intens;
}

const ParameterList& AmpSumIntensity::intensityInterference(dataPoint& point,
		resonanceItr A, resonanceItr B)
{
	double intens = (
			(*A)->Evaluate(point)*std::conj((*B)->Evaluate(point))
	).real();
	if( A != B ) intens = 2*intens;
	result.SetParameterValue(0,intens);

	return result;
}

const double AmpSumIntensity::interferenceIntegral(
		resonanceItr A, resonanceItr B)
{
	tmpA = A; tmpB = B;

	/* Integration functionality was tested with a model with only one
	 * normalized amplitude. The integration result is equal to the
	 * amplitude coefficient^2.
	 */
	size_t dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	auto m13sq_limit = kin->GetMinMax(1);
	auto m23sq_limit = kin->GetMinMax(0);
	double xLimit_low[2] = {m13sq_limit.first,m23sq_limit.first};
	double xLimit_high[2] = {m13sq_limit.second,m23sq_limit.second};
	//	double xLimit_low[2] = {0,0};
	//	double xLimit_high[2] = {10,10};
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function G = {&interferenceIntegralWrapper,dim, const_cast<AmpSumIntensity*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&G, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::interferenceIntegrate() | "
			"Interference term of "<<(*tmpA)->GetName()<<" and "<<(*tmpB)->GetName()
			<<" : "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

double AmpSumIntensity::evaluateEff(double x[], size_t dim)
{
	/* Calculation amplitude normalization (including efficiency) */
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//set data point: we assume that x[0]=m13 and x[1]=m23
	dataPoint point;
	try{
		kin->FillDataPoint(1,0,x[0],x[1],point);
	} catch (BeyondPhsp& ex){
		return 0;
	}
	ParameterList res = intensity(point);
	double intens = *res.GetDoubleParameter(0);
	return intens;
}

double evalWrapperAmpSumIntensityEff(double* x, size_t dim, void* param)
{
	/* We need a wrapper here because intensity() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	return static_cast<AmpSumIntensity*>(param)->evaluateEff(x,dim);
}

const double AmpSumIntensity::normalization(ParameterList& par)
{
	setParameterList(par);
	return normalization();
}
const double AmpSumIntensity::normalization()
{
	/* Integration functionality was tested with a model with only one normalized amplitude.
	 * The integration result is equal to the amplitude coefficient^2.
	 */
	size_t dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	auto m13sq_limit = kin->GetMinMax(1);
	auto m23sq_limit = kin->GetMinMax(0);
	double xLimit_low[2] = {m13sq_limit.first,m23sq_limit.first};
	double xLimit_high[2] = {m13sq_limit.second,m23sq_limit.second};
	//	double xLimit_low[2] = {0,0};
	//	double xLimit_high[2] = {10,10};
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function G = {&evalWrapperAmpSumIntensityEff,dim, const_cast<AmpSumIntensity*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&G, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	//	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity::normalization() Integration result for amplitude sum: "
	//			<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

const double AmpSumIntensity::sliceIntensity(dataPoint& dataP,
		ParameterList& par,std::complex<double>* reso, unsigned int nResos)
{
	setParameterList(par);

	double AMPpdf=0;
	//TODO: implement slice fit
	//	if(Kinematics::instance()->isWithinPhsp(dataP)) AMPpdf = totAmp.evaluateSlice(dataP, reso, nResos,5);
	if(AMPpdf!=AMPpdf){
		BOOST_LOG_TRIVIAL(error)<<"Error AmpSumIntensity: Intensity is not a number!!";
		AMPpdf = 0;
	}
	double eff=eff_->evaluate(dataP);
	return AMPpdf*eff;
}

const ParameterList& AmpSumIntensity::intensity(std::vector<double> point,
		ParameterList& par)
{
	setParameterList(par);
	dataPoint dataP(point);
	return intensity(dataP);
}

const ParameterList& AmpSumIntensity::intensity(dataPoint& point, ParameterList& par)
{
	setParameterList(par);
	return intensity(point);
}

const ParameterList& AmpSumIntensity::intensityNoEff(dataPoint& point)
{
	std::complex<double> AMPpdf(0,0);
	if(Kinematics::instance()->isWithinPhsp(point)) {
		auto it = GetResonanceItrFirst();
		for( ; it != GetResonanceItrLast(); ++it){
			AMPpdf += (*it)->Evaluate(point);
		}
	}
	result.SetParameterValue(0,std::norm(AMPpdf));
	return result;
}

const ParameterList& AmpSumIntensity::intensity(dataPoint& point)
{
	intensityNoEff(point);
	double ampNoEff = result.GetParameterValue(0);
	double eff=eff_->evaluate(point);
	result.SetParameterValue(0,ampNoEff*eff);
	return result;
}

void AmpSumIntensity::setParameterList(ParameterList& par)
{
	//parameters varied by Minimization algorithm
	if(par.GetNDouble()!=params.GetNDouble())
		throw std::runtime_error("AmpSumIntensity::setParameterList(): size of parameter lists don't match");
	//Should we compared the parameter names? String comparison is slow
	for(unsigned int i=0; i<params.GetNDouble(); i++)
		params.GetDoubleParameter(i)->UpdateParameter(par.GetDoubleParameter(i));
	return;
}

bool AmpSumIntensity::copyParameterList(ParameterList& outPar)
{
	outPar = ParameterList(params);
	return true;
}

void AmpSumIntensity::printAmps()
{
	std::stringstream outStr;
	outStr<<"AmpSumIntensity: Printing amplitudes with current(!) set of parameters:\n";
	unsigned int n=0;
	for(unsigned int i=0; i<params.GetNDouble(); i++){
		std::shared_ptr<DoubleParameter> p = params.GetDoubleParameter(i);
		std::string tmp = p->GetName();
		std::size_t found = tmp.find("mag");
		if(found!=std::string::npos){
			outStr<<"-------- "<<GetNameOfResonance(n)<<" ---------\n";
			n++;
		}
		outStr<<p->GetName()<<" = "<<p->GetValue();
		if(p->HasError())
			outStr<<"+-"<<p->GetError();
		if(p->HasBounds())
			outStr<<" ["<<p->GetMinValue()<<";"<<p->GetMaxValue()<<"]";
		if(p->IsFixed())
			outStr<<" FIXED";
		outStr<<"\n";
	}

	BOOST_LOG_TRIVIAL(info)<<outStr.str();
	return;
}

void AmpSumIntensity::printFractions()
{
	std::stringstream outStr;
	outStr<<"Fit fractions for all amplitudes: \n";
	double sumFrac=0;
	auto it = GetResonanceItrFirst();
	for( ; it != GetResonanceItrLast(); ++it){
		double frac = (*it)->GetMagnitude()/integral();
		sumFrac+=frac;
		outStr<<std::setw(10)<<(*it)->GetName()<<":    "<<frac<<"\n";
	}

	outStr<<std::setw(10)<<" "<<"    ==========\n";
	outStr<<std::setw(10)<<" "<<"     "<<sumFrac;
	BOOST_LOG_TRIVIAL(info)<<outStr.str();
	return;
}

double AmpSumIntensity::getIntValue(std::string var1, double min1, double max1,
		std::string var2, double min2, double max2)
{
	/*
	 * Integrates in \var1 from \min1 to \max1 and in \var2 from \min2 to \max2.
	 * Is intended to be used for calculation of bin content.
	 */
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double _min1 = min1;
	double _min2 = min2;
	double _max1 = max1;
	double _max2 = max2;
	//if(_min1==0) _min1 = kin->GetMin(var1);
	//if(_max1==0) _min1 = kin->GetMax(var1);
	auto limit2 = kin->GetMinMax(var2);
	if(_min2==0) _min2 = limit2.first;
	if(_max2==0) _max2 = limit2.second;
	unsigned int dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	double xLimit_low[2];
	double xLimit_high[2];

	if(var1 == "m13sq" && var2 == "m23sq"){
		xLimit_low[0] = _min1;
		xLimit_low[1] = _min2;
		xLimit_high[0] = _max1;
		xLimit_high[1] = _max2;
	}else if(var1 == "m23sq" && var2 == "m13sq"){
		xLimit_low[0] = _min2;
		xLimit_low[1] = _min1;
		xLimit_high[0] = _max2;
		xLimit_high[1] = _max1;
	} else {
		BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::getIntValue() wrong variables specified!";
		return -999;
	}

	size_t calls = 5000;
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function G = {&evalWrapperAmpSumIntensity,dim, const_cast<AmpSumIntensity*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&G, xLimit_low, xLimit_high, 2, calls, r,s,&res, &err);
	gsl_monte_vegas_free(s);

	return res;
}

int AmpSumIntensity::GetIdOfResonance(std::string name)
{
	for(unsigned int i=0; i< resoList.size(); i++)
		if(resoList.at(i)->GetName()==name) return i;
	return -999;
}

std::string AmpSumIntensity::GetNameOfResonance(unsigned int id)
{
	if(id < 0 || id > resoList.size() )
		throw std::runtime_error("AmpSumIntensity::getNameOfResonance() | "
				"Invalid resonance ID! Resonance not found?");
	return resoList.at(id)->GetName();
}

std::shared_ptr<Resonance> AmpSumIntensity::GetResonance(std::string name)
{
	int id = GetIdOfResonance(name);
	return GetResonance(id);
}

std::shared_ptr<Resonance> AmpSumIntensity::GetResonance(unsigned int id)
{
	if(id < 0 || id > resoList.size() )
		throw std::runtime_error("AmpSumIntensity::getResonance() | "
				"Invalid resonance ID! Resonance not found?");
	return resoList.at(id);
}

double AmpSumIntensity::averageWidth()
{
	double avWidth = 0;
	double sum = 0;
	for(int i=0; i<resoList.size(); i++){
		avWidth += std::norm(resoList.at(i)->GetMagnitude())*resoList.at(i)->GetWidth();
		sum += std::norm(resoList.at(i)->GetMagnitude());
	}
	avWidth /= sum;
	return avWidth;
}
