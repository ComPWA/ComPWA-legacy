#include "DataReader/Data.hpp"

void rndReduceSet(unsigned int size, std::shared_ptr<Generator> gen,
		 Data* in1, Data* out1, Data*  in2, Data*  out2){
	if(!in1 )
		throw std::runtime_error("rndSubSet() | No input data set!");
	if(!out1 )
		throw std::runtime_error("rndSubSet() | No output data set!");
	if( out1->getNEvents() )
		throw std::runtime_error("rndSubSet() | First output sample not empty!");
	if(in2){
		std::cout<<in1->getNEvents() <<" "<< in2->getNEvents()<<std::endl;
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
		dataPoint point(in1->getEvent(i));
		if(!Kinematics::instance()->isWithinPhsp(point)) newSize--;
	}
	double threshold = (double)size/newSize; //calculate threshold
	for(unsigned int i=0; i<totalSize; i++){
		dataPoint point(in1->getEvent(i)); //use first sample for hit&miss
		if(!Kinematics::instance()->isWithinPhsp(point)) continue;
		if( gen->getUniform() < threshold ){
			out1->pushEvent(in1->getEvent(i));
			//write second sample if event from first sample were accepted
			if(in2) out2->pushEvent(in2->getEvent(i));
		}
	}
	std::cout<<out1->getNEvents()<<std::endl;
	std::cout<<"asdfasdfasfasdfsdfasdfsad"<<std::endl;
	return;
}

std::shared_ptr<Data> Data::rndSubSet(unsigned int size, std::shared_ptr<Generator> gen){
	std::shared_ptr<Data> out(this->EmptyClone());
	rndReduceSet(size, gen,this, out.get());
	return out;
}
