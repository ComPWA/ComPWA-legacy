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

#include <omp.h>
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

#include "Core/RunManager.hpp"
using namespace boost::log;

RunManager::RunManager() : success_(false), validAmplitude(0), validOptimizer(0),
		validBackground(0), validGenerator(0), size_(0), bkgSize_(0)
{
}
RunManager::RunManager(std::shared_ptr<Data> inD, std::shared_ptr<Amplitude> inP,
		std::shared_ptr<Optimizer> inO) : pData_(inD), pPhys_(inP), pOpti_(inO), success_(false),
  validAmplitude(0), validOptimizer(0), validBackground(0), bkgSize_(0)
{
	setSize(inD->getNEvents());
	if(inD && inP && inO){
		validAmplitude=1;
		validOptimizer=1;
	}
}
RunManager::RunManager( unsigned int size, std::shared_ptr<Amplitude> inP,
		std::shared_ptr<Generator> gen) : gen_(gen), size_(size), pPhys_(inP), success_(false),
  validAmplitude(0), validOptimizer(0), validBackground(0), bkgSize_(0)
{
	if(inP ){
		validAmplitude=1;
	}
}

RunManager::~RunManager(){
	BOOST_LOG_TRIVIAL(info) << "~RunManager: Last seed: "<<gen_->getSeed();
}

std::shared_ptr<FitResult> RunManager::startFit(ParameterList& inPar){
	BOOST_LOG_TRIVIAL(info) << "RunManager: Starting fit.";
	BOOST_LOG_TRIVIAL(info) << "RunManager: Input data contains "<<pData_->getNEvents()<<" events.";
	std::shared_ptr<FitResult> result = pOpti_->exec(inPar);
	success_ = true;
	BOOST_LOG_TRIVIAL(info) << "RunManager: fit end. Result ="<<result->getResult()<<".";

	return result;
}
bool RunManager::generate( unsigned int number ) {
	if(number>0) setSize(number);

	if( !(pData_ && validAmplitude && validGenerator) )
		throw std::runtime_error("RunManager: generate() requirements not fulfilled");
	if(pData_->getNEvents()>0)
		throw std::runtime_error("RunManager: generate() dataset not empty! abort!");

	//Determing an estimate on the maximum of the physics amplitude using 100k events.
	double genMaxVal=1.2*pPhys_->getMaxVal(gen_);

	unsigned int totalCalls=0;
	BOOST_LOG_TRIVIAL(info) << "Generating MC: ["<<size_<<" events] ";
	BOOST_LOG_TRIVIAL(info) << "== Using "<<genMaxVal<< " as maximum value for random number generation!";

	unsigned int startTime = clock();

	Generator* genNew = (&(*gen_))->Clone();//copy generator for every thread
	double AMPpdf;

	unsigned int cnt=0;
	for(unsigned int i=0;i<size_;i++){
		if(i>0) i--;
		if(i%1000==0 && cnt!=i) {
			cnt=i;
			BOOST_LOG_TRIVIAL(debug) << "Current event: "<<i;
		}
		Event tmp;
		genNew->generate(tmp);
		totalCalls++;
		double weight = tmp.getWeight();
		/* reset weights: the weights are taken into account by hit and miss. The resulting
		 * sample is therefore unweighted */
		tmp.setWeight(1.);//reset weight
		tmp.setEfficiency(1.);//reset weight
		dataPoint point(tmp);
		double ampRnd = genNew->getUniform()*genMaxVal;
		ParameterList list;
		list = pPhys_->intensity(point);//unfortunatly not thread safe
		AMPpdf = *list.GetDoubleParameter(0);
		if( genMaxVal< (AMPpdf*weight))
			throw std::runtime_error("RunManager: error in HitMiss procedure. "
					"Maximum value of random number generation smaller then amplitude maximum!");
		if( ampRnd > (weight*AMPpdf) ) continue;
		i++;
		pData_->pushEvent(tmp);//Unfortunately not thread safe
	}
//	BOOST_LOG_TRIVIAL(info) << "Efficiency of toy MC generation: "<<(double)size_/totalCalls;
//	BOOST_LOG_TRIVIAL(info) << "RunManager: generate time="<<(clock()-startTime)/CLOCKS_PER_SEC/60<<"min.";

	return true;
};
bool RunManager::generateBkg( unsigned int number ) {
	if(number>0) setBkgSize(number);
	if(!bkgSize_) return true;

	if( !(ampBkg_ && validBackground && validGenerator) )
		throw std::runtime_error("RunManager: generateBkg() requirements not fulfilled");
	if(sampleBkg_->getNEvents()>0)
		throw std::runtime_error("RunManager: generateBkg() dataset not empty! abort!");
	//Determing an estimate on the maximum of the physics amplitude using 100k events.
	double genMaxVal=1.2*ampBkg_->getMaxVal(gen_);

	unsigned int totalCalls=0;
	BOOST_LOG_TRIVIAL(info) << "Generating background MC: ["<<bkgSize_<<" events] ";

	unsigned int startTime = clock();

	Generator* genNew = (&(*gen_))->Clone();//copy generator for every thread
	double AMPpdf;

	unsigned int cnt=0;
	for(unsigned int i=0;i<bkgSize_;i++){
		if(i>0) i--;
		if(i%1000==0 && cnt!=i) {
			cnt=i;
			BOOST_LOG_TRIVIAL(debug) << "Current event: "<<i;
		}
		Event tmp;
		genNew->generate(tmp);
		totalCalls++;
		double weight = tmp.getWeight();
		/* reset weights: the weights are taken into account by hit and miss. The resulting
		 * sample is therefore unweighted */
		tmp.setWeight(1.);//reset weight
		tmp.setEfficiency(1.);//reset weight
		dataPoint point(tmp);
		double ampRnd = genNew->getUniform()*genMaxVal;
		ParameterList list;
		list = ampBkg_->intensity(point);//unfortunatly not thread safe
		AMPpdf = *list.GetDoubleParameter(0);
		if( genMaxVal< (AMPpdf*weight))
			throw std::runtime_error("RunManager: error in HitMiss procedure. "
					"Maximum value of random number generation smaller then amplitude maximum!");
		if( ampRnd > (weight*AMPpdf) ) continue;
		i++;
		sampleBkg_->pushEvent(tmp);//unfortunatly not thread safe
	}
//	BOOST_LOG_TRIVIAL(info) << "RunManager: generate background time="<<(clock()-startTime)/CLOCKS_PER_SEC/60<<"min.";

	return true;
};

bool RunManager::generatePhsp( unsigned int number ) {
	if( !(validPhsp==1) )
		return false;
	unsigned int phspSize = 10*size_;
	if(number>0) phspSize = number;

	BOOST_LOG_TRIVIAL(info) << "Generating phase-space MC: ["<<phspSize<<" events] ";

	Generator* genNew = (&(*gen_))->Clone();//copy generator for every thread

	for(unsigned int i=0;i<phspSize;i++){
		if(i>0) i--;
		Event tmp;
		genNew->generate(tmp);
		double ampRnd = genNew->getUniform();
		if( ampRnd > tmp.getWeight() ) continue;
		/* reset weights: the weights are taken into account by hit and miss on the weights.
		 * The resulting sample is therefore unweighted */
		tmp.setWeight(1.);//reset weight
		tmp.setEfficiency(1.);
		i++;
		phspSample_->pushEvent(tmp);//unfortunatly not thread safe
	}
	return true;
};
