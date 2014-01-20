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
//#include <boost/log/core.hpp>
//#include <boost/log/expressions.hpp>
#include <boost/progress.hpp>

#include "DataReader/Data.hpp"
#include "Estimator/Estimator.hpp"
#include "Physics/Amplitude.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Event.hpp"
#include "Core/Generator.hpp"

#include "Core/RunManager.hpp"
//#define BOOST_LOG_DYN_LINK 1
using namespace boost::log;

RunManager::RunManager(std::shared_ptr<Data> inD, std::shared_ptr<Amplitude> inP,
		std::shared_ptr<Optimizer> inO)
: pData_(inD), pPhys_(inP), pOpti_(inO), success_(false),
  validSize(0),validAmplitude(0),validData(0),validOptimizer(0)
{
	setSize(inD->getNEvents());
	if(inD && inP && inO){
		validAmplitude=1;
		validData=1;
		validOptimizer=1;
	}
}
RunManager::RunManager( unsigned int size, std::shared_ptr<Amplitude> inP,
		std::shared_ptr<Generator> gen)
: gen_(gen), size_(size), pPhys_(inP), success_(false),
  validSize(0),validAmplitude(0),validData(0),validOptimizer(0)
{
	if(inP ){
		validAmplitude=1;
	}
	validSize=1;
}

RunManager::~RunManager(){
	/* nothing */
}

std::shared_ptr<FitResult> RunManager::startFit(ParameterList& inPar){
	if( !(validAmplitude==1 && validData==1 && validOptimizer==1) ){
		BOOST_LOG_TRIVIAL(error) << "RunManager: startFit() requirements not fulfilled!";
		return std::shared_ptr<FitResult>();
	}
	BOOST_LOG_TRIVIAL(info) << "RunManager: Starting fit.";
	BOOST_LOG_TRIVIAL(info) << "RunManager: Input data contains "<<pData_->getNEvents()<<" events.";
	//	BOOST_LOG_TRIVIAL(info) << "Initial LH="<<esti->controlParameter(inPar);
	unsigned int startTime = clock();
	std::shared_ptr<FitResult> result = pOpti_->exec(inPar);
	success_ = true;
	BOOST_LOG_TRIVIAL(info) << "RunManager: fit end. Result ="<<result->getResult()<<".";
	BOOST_LOG_TRIVIAL(info) << "RunManager: fit time="<<(clock()-startTime)/CLOCKS_PER_SEC/60<<"min.";

	return result;
}
bool RunManager::generate( unsigned int number ) {
	if(number>0) setSize(number);
	if( !(validData && validAmplitude && validSize && validGenerator) ){
		BOOST_LOG_TRIVIAL(error)<<"RunManager: generate() requirements not fulfilled";
		return false;
	}

	if(pData_->getNEvents()>0) {
		BOOST_LOG_TRIVIAL(error)<<"RunManager: generate() dataset not empty! abort!";
		return false;
	}
	ParameterList minPar;
	pPhys_->fillStartParVec(minPar);

	//Determing an estimate on the maximum of the physics amplitude using 20k events.
	//	double genMaxVal=2*pPhys_->getMaxVal();
	double genMaxVal=pPhys_->getMaxVal();

	double maxTest=0;
	unsigned int totalCalls=0;
	BOOST_LOG_TRIVIAL(info) << "== Using "<<genMaxVal<< " as maximum value for random number generation!";
	BOOST_LOG_TRIVIAL(info) << "Generating MC: ["<<size_<<" events] ";

	unsigned int startTime = clock();
	boost::progress_display progressBar(size_); //boost progress bar (thread-safe)
#pragma omp parallel firstprivate(genMaxVal) shared(maxTest,totalCalls)
	{
		unsigned int threadId = omp_get_thread_num();
		unsigned int numThreads = omp_get_num_threads();

		Generator* genNew = (&(*gen_))->Clone();//copy generator for every thread
		//		genNew->setSeed(std::clock()+threadId);//setting the seed here makes not sense in cast that TGenPhaseSpace is used, because it uses gRandom
		//initialize random number generator
		boost::minstd_rand rndGen2(std::clock()+threadId);//TODO: is this seed thread safe?
		boost::uniform_real<> uni_dist2(0,1);
		boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uni2(rndGen2, uni_dist2);
		double AMPpdf;

#pragma omp for
		for(unsigned int i=0;i<size_;i++){
			if(i>0) i--;
			//	while( i<size_){ //while loops are not supported by openMP
			Event tmp;
			genNew->generate(tmp);
			totalCalls++;
			double weight = tmp.getWeight();
			tmp.setWeight(1);//reset weight
//			Particle part1 = tmp.getParticle(0);
//			Particle part2 = tmp.getParticle(1);
//			Particle part3 = tmp.getParticle(2);
//			double m23sq = Particle::invariantMass(part2,part3);
//			double m13sq = Particle::invariantMass(part1,part3);
//			std::vector<double> x;
//			x.push_back(m23sq);
//			x.push_back(m13sq);
//			std::cout<<point.getVal("m23sq")<<" "<<x[0]<< " "<<point.getVal("m13sq")<<" "<<x[1]<<std::endl;
			dataPoint point(tmp);

			double ampRnd = uni2()*genMaxVal;
			ParameterList list;
#pragma omp critical
			{
				list = pPhys_->intensity(point,minPar);//unfortunatly not thread safe
				AMPpdf = *list.GetDoubleParameter(0);
				if( maxTest < (AMPpdf*weight)) maxTest = weight*AMPpdf;
			}
			if( ampRnd > (weight*AMPpdf) ) continue;
			i++;
#pragma omp critical
			{
				pData_->pushEvent(tmp);//unfortunatly not thread safe
			}
			++progressBar;//progress bar
		}
	}
	BOOST_LOG_TRIVIAL(info) << "Max value used in generation: "<<maxTest;
	BOOST_LOG_TRIVIAL(info) << "Efficiency of toy MC generation: "<<(double)size_/totalCalls;
	BOOST_LOG_TRIVIAL(info) << "RunManager: generate time="<<(clock()-startTime)/CLOCKS_PER_SEC/60<<"min.";

	if( maxTest > (double) (0.9*genMaxVal) ) {
		BOOST_LOG_TRIVIAL(error)<<"==========ATTENTION===========";
		BOOST_LOG_TRIVIAL(error)<<"== Max value of function is "<<maxTest;
		BOOST_LOG_TRIVIAL(error)<<"== This is close or above to maximum value of rnd. number generation: "<<genMaxVal;
		BOOST_LOG_TRIVIAL(error)<<"== Choose higher max value!";
		BOOST_LOG_TRIVIAL(error)<<"==========ATTENTION===========";
		return false;
	}

	return true;
};
bool RunManager::generatePhsp( unsigned int number ) {
	if( !(validPhsp==1 && validSize==1) )
		return false;
	unsigned int phspSize = 10*size_;
	if(number>0) phspSize = number;

	BOOST_LOG_TRIVIAL(info) << "Generating phase-space MC: ["<<phspSize<<" events] ";

//	boost::progress_display progressBar(phspSize); //boost progress bar (thread-safe)
#pragma omp parallel
	{
		Generator* genNew = (&(*gen_))->Clone();//copy generator for every thread
		//		genNew->setSeed(std::clock()+threadId);//setting the seed here makes not sense in cast that TGenPhaseSpace is used, because it uses gRandom

#pragma omp for
		for(unsigned int i=0;i<phspSize;i++){
			Event tmp;
			genNew->generate(tmp);
#pragma omp critical
			{
				phspSample_->pushEvent(tmp);//unfortunatly not thread safe
			}
//			++progressBar;//progress bar
		}
	}
	return true;
};
