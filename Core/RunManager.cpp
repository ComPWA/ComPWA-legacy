//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding functionality to generate set of events
//-------------------------------------------------------------------------------
#include <memory>
#include <ctime>

#include <time.h>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/progress.hpp>

#include "DataReader/Data.hpp"
#include "Estimator/Estimator.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Event.hpp"
#include "Core/Generator.hpp"
#include "Core/ProgressBar.hpp"

#include "Core/RunManager.hpp"
using namespace boost::log;

RunManager::RunManager()
{
}
RunManager::RunManager(std::shared_ptr<Data> inD, std::shared_ptr<Amplitude> inP,
		std::shared_ptr<Optimizer> inO) : sampleData_(inD), amp_(inP), opti_(inO)
{
}
RunManager::RunManager( unsigned int size, std::shared_ptr<Amplitude> inP,
		std::shared_ptr<Generator> gen) : gen_(gen),  amp_(inP)
{
}

RunManager::~RunManager(){
	BOOST_LOG_TRIVIAL(info) << "~RunManager: Last seed: "<<gen_->getSeed();
}

std::shared_ptr<FitResult> RunManager::startFit(ParameterList& inPar){
	BOOST_LOG_TRIVIAL(info) << "RunManager: Starting fit.";
	BOOST_LOG_TRIVIAL(info) << "RunManager: Input data contains "<<sampleData_->getNEvents()<<" events.";
	std::shared_ptr<FitResult> result = opti_->exec(inPar);
	BOOST_LOG_TRIVIAL(info) << "RunManager: fit end. Result ="<<result->getResult()<<".";

	return result;
}

void RunManager::setPhspSample( std::shared_ptr<Data> phsp, std::shared_ptr<Data> truePhsp){
	if(truePhsp && truePhsp->getNEvents() != phsp->getNEvents())
		throw std::runtime_error("RunManager::setPhspSample() | "
				"Reconstructed sample and true sample have not the same size!");
	samplePhsp_ = phsp;
	sampleTruePhsp_ = truePhsp;
}
void RunManager::setTruePhspSample( std::shared_ptr<Data> truePhsp){
	if(truePhsp && samplePhsp_ && truePhsp->getNEvents() != samplePhsp_->getNEvents())
		throw std::runtime_error("RunManager::setPhspSample() | "
				"Reconstructed sample and true sample have not the same size!");
	sampleTruePhsp_ = truePhsp;
};

bool RunManager::generate( int number ) {
	BOOST_LOG_TRIVIAL(info) << "RunManager::generate() | generating "<<number<<" signal events!";
	return gen( number, gen_, amp_, sampleData_, samplePhsp_, sampleTruePhsp_ );
}

bool RunManager::generateBkg( int number ) {
	BOOST_LOG_TRIVIAL(info) << "RunManager::generateBkg() | generating "<<number<<" background events!";
	return gen( number, gen_, ampBkg_, sampleBkg_, samplePhsp_, sampleTruePhsp_);
}

bool RunManager::gen( int number, std::shared_ptr<Generator> gen, std::shared_ptr<Amplitude> amp,
		std::shared_ptr<Data> data, std::shared_ptr<Data> phsp, std::shared_ptr<Data> phspTrue){

	if(number == 0) return 0;

	//Doing some checks
	if(number < 0 && !phsp)
		throw std::runtime_error("RunManager: gen() negative number of events: "
				+std::to_string((long double)number)+". And no phsp sample given!");
	if( !amp )
		throw std::runtime_error("RunManager: gen() amplitude not valid");
	if( !gen )
		throw std::runtime_error("RunManager: gen() generator not valid");
	if( !data )
		throw std::runtime_error("RunManager: gen() sample not valid");
	if( data->getNEvents() > 0 )
		throw std::runtime_error("RunManager: gen() sample not empty!");
	if( phspTrue && !phsp )
		throw std::runtime_error("RunManager: gen() We have a sample of true phsp events, "
				"but no phsp sample!");
	if( phspTrue && phspTrue->getNEvents() != phsp->getNEvents() )
		throw std::runtime_error("RunManager: gen() We have a sample of true phsp events, "
				"but the sample size doesn't match that one of the phsp sample!");

	double maxSampleWeight = 1.0;
	if(phsp) maxSampleWeight = phsp->getMaxWeight();
	if(phspTrue && phspTrue->getMaxWeight()>maxSampleWeight) maxSampleWeight = phspTrue->getMaxWeight();

	/* Maximum value for random number generation. We introduce an arbitrary factor of 1.2 to make sure
	 * that the maximum value is never reached.
	 */
	double generationMaxValue = 1.2*amp->GetMaxVal(gen)*maxSampleWeight;

	unsigned int totalCalls=0;
	unsigned int startTime = clock();
	double AMPpdf;

	unsigned int limit;
	unsigned int acceptedEvents=0;
	if(phsp){
		limit = phsp->getNEvents();
	} else
		limit = 100000000; //set large limit, should never be reached

	Event evt; //event that we fill into generated sample
	Event evtTrue; //event that is used to evalutate amplitude
	progressBar bar(number);
	if(number<=0) bar = progressBar(limit);
	for(unsigned int i=0;i<limit;i++){
		if(phsp && phspTrue) { //phsp and true sample is set
			evtTrue = phspTrue->getEvent(i);
			evt = phsp->getEvent(i);
		} else if(phsp) {//phsp sample is set
			evt = phsp->getEvent(i);
			evtTrue = evt;
		} else {//otherwise generate event
			gen->generate(evt);
			evtTrue = evt;
		}
		if(number<=0) bar.nextEvent();

		//use reconstructed position for weights
		double weight = evt.getWeight();

		//use true position for amplitude value
		dataPoint point;
		try{
			point = dataPoint(evtTrue);
		} catch (BeyondPhsp& ex){ //event outside phase, remove
			continue;
		}
		//		if(!Kinematics::instance()->isWithinPhsp(point)) continue;

		totalCalls++;
		double ampRnd = gen->getUniform()*generationMaxValue;
		ParameterList list;
		list = amp->intensity(point);//unfortunatly not thread safe
		AMPpdf = *list.GetDoubleParameter(0);
		if( generationMaxValue< (AMPpdf*weight))
			throw std::runtime_error("RunManager: error in HitMiss procedure. "
					"Maximum value of random number generation smaller then amplitude maximum!");
		if( ampRnd > (weight*AMPpdf) ) continue;

		//Fill event to sample
		/* reset weights: the weights are taken into account by hit and miss. The resulting
		 * sample is therefore unweighted */
		evt.setWeight(1.);//reset weight
		evt.setEfficiency(1.);//reset weight
		data->pushEvent(evt);//unfortunatly not thread safe

		//some statistics
		acceptedEvents++;
		if(number>0) bar.nextEvent();
		if(acceptedEvents>=number) i=limit; //break if we have a sufficienct number of events
	}
	if(data->getNEvents()<number)
		BOOST_LOG_TRIVIAL(error) << "RunManager::gen() not able to generate "<<number<<" events. "
				"Phsp sample too small. Current size of sample is now "<<data->getNEvents();
	BOOST_LOG_TRIVIAL(info) << "Efficiency of toy MC generation: "<<(double)data->getNEvents()/totalCalls;

	return true;
}

bool RunManager::generatePhsp( int number ) {
	if(number==0) return 0;
	if(!samplePhsp_)
		throw std::runtime_error("RunManager: generatePhsp() not phsp sample set");
	if(samplePhsp_->getNEvents()>0)
		throw std::runtime_error("RunManager: generatePhsp() dataset not empty! abort!");

	BOOST_LOG_TRIVIAL(info) << "Generating phase-space MC: ["<<number<<" events] ";

	progressBar bar(number);
	for(unsigned int i=0;i<number;i++){
		if(i>0) i--;
		Event tmp;
		gen_->generate(tmp);
		double ampRnd = gen_->getUniform();
		if( ampRnd > tmp.getWeight() ) continue;
		/* reset weights: the weights are taken into account by hit and miss on the weights.
		 * The resulting sample is therefore unweighted */
		tmp.setWeight(1.);//reset weight
		tmp.setEfficiency(1.);
		i++;
		samplePhsp_->pushEvent(tmp);//unfortunatly not thread safe
		bar.nextEvent();
	}
	return true;
}
