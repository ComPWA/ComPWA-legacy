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

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

#include "DataReader/Data.hpp"
#include "Estimator/Estimator.hpp"
#include "Physics/Amplitude.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Event.hpp"
#include "Core/Generator.hpp"

#include "Core/RunManager.hpp"

RunManager::RunManager(std::shared_ptr<Data> inD, std::shared_ptr<Amplitude> inP,
		std::shared_ptr<Optimizer> inO, std::shared_ptr<Efficiency> eff)
: eff_(eff), pData_(inD), pPhys_(inP), pOpti_(inO), success_(false),
  validSize(0),validAmplitude(0),validData(0),validOptimizer(0)
{
	if(eff && inD && inP && inO){
		validEfficiency=1;
		validAmplitude=1;
		validData=1;
		validOptimizer=1;
	}
}
RunManager::RunManager( unsigned int size, std::shared_ptr<Amplitude> inP,
		std::shared_ptr<Efficiency> eff, std::shared_ptr<Generator> gen)
: gen_(gen), eff_(eff), size_(size), pPhys_(inP), success_(false),
  validSize(0),validAmplitude(0),validData(0),validOptimizer(0)
{
	if(inP && eff){
		validEfficiency=1;
		validAmplitude=1;
	}
	validSize=1;
}

RunManager::~RunManager(){
	/* nothing */
}

bool RunManager::startFit(ParameterList& inPar){
	if( !(validEfficiency==1 && validAmplitude==1 && validData==1 && validOptimizer==1) )
		return false;

	pOpti_->exec(inPar);
	success_ = true;

	return success_;
}
bool RunManager::generate( unsigned int number ) {
	if( !(validData==1 && validEfficiency==1 && validAmplitude==1 && validSize==1) )
		return false;

	//	if(pData_->getNEvents()>0){
	//What do we do if dataset is not empty?
	//	}
	ParameterList minPar;
	pPhys_->fillStartParVec(minPar);

	//initialize random number generator
	boost::minstd_rand rndGen(std::clock());//TODO: is this seed thread safe?
	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uni(rndGen, uni_dist);

	//Determing an estimate on the maximum of the physics amplitude using 10k events.
	double genMaxVal=0;
	for(unsigned int i=0; i<20000;i++){
		Event tmp;
		gen_->generate(tmp);
		double weight = tmp.getWeight();
		Particle part1 = tmp.getParticle(0);
		Particle part2 = tmp.getParticle(1);
		Particle part3 = tmp.getParticle(2);
		double m23sq = Particle::invariantMass(part2,part3);
		double m13sq = Particle::invariantMass(part1,part3);
		double m12sq = Particle::invariantMass(part1,part2);
		std::vector<double> x;
		x.push_back(m23sq);
		x.push_back(m13sq);
		x.push_back(m12sq);
		//		double eff  = eff_->evaluate(x);
		ParameterList list = pPhys_->intensity(x,minPar);
		double AMPpdf = *list.GetDoubleParameter(0);
		if(genMaxVal<(weight*AMPpdf)) genMaxVal= weight*AMPpdf;
	}
	genMaxVal=2*genMaxVal;//conservative choice

	unsigned int i=0;
	double maxTest=0;
	int scale = (int) size_/10;
	std::cout<<"== Using "<<genMaxVal<< " as maximum value for random number generation!"<<std::endl;
	std::cout << "Generating MC: ["<<size_<<" events] ";
	while( i<size_){
		Event tmp;
		gen_->generate(tmp);
		double weight = tmp.getWeight();
		Particle part1 = tmp.getParticle(0);
		Particle part2 = tmp.getParticle(1);
		Particle part3 = tmp.getParticle(2);
		double m23sq = Particle::invariantMass(part2,part3);
		double m13sq = Particle::invariantMass(part1,part3);
		double m12sq = Particle::invariantMass(part1,part2);
		std::vector<double> x;
		x.push_back(m23sq);
		x.push_back(m13sq);
		x.push_back(m12sq);
		double eff  = eff_->evaluate(x);
		//Efficiency is saved to event. Weightng is done when parameters are plotted.
		tmp.setWeight(eff);

		double ampRnd = uni()*genMaxVal;
		ParameterList list = pPhys_->intensity(x,minPar);
		double AMPpdf = *list.GetDoubleParameter(0);

		//		if( maxTest < (AMPpdf*weight*eff)) maxTest = AMPpdf;
		//		if( ampRnd>(weight*AMPpdf*eff) ) continue;
		if( maxTest < (AMPpdf*weight)) maxTest = weight*AMPpdf;
		if( ampRnd > (weight*AMPpdf) ) continue;
		pData_->pushEvent(tmp);
		i++;

		//progress bar
		if( (i % scale) == 0) { std::cout<<(i/scale)*10<<"%..."<<std::flush; }
	}
	std::cout<<std::endl;

	if( maxTest > (int) (0.9*genMaxVal) ) {
		std::cout<<"==========ATTENTION==========="<<std::endl;
		std::cout<<"== Max value of function is "<<maxTest<<std::endl;
		std::cout<<"== This is close or above to maximum value of rnd. number generation: "<<genMaxVal<<std::endl;
		std::cout<<"== Choose higher max value!"<<std::endl;
		std::cout<<"==========ATTENTION==========="<<std::endl;
		return false;
	}

	return true;
};
