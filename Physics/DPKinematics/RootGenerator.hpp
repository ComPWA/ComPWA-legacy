//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

#ifndef ROOTGENERATOR_HPP_
#define ROOTGENERATOR_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

#include "TLorentzVector.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"

#include "Core/Generator.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

class RootGenerator : public Generator {
protected:
	TGenPhaseSpace event;
	unsigned int nPart;
	Double_t* masses;
public:
	RootGenerator(int seed=-1) ;
	~RootGenerator(){};

	virtual RootGenerator* Clone();
	virtual void generate(Event& evt);
	virtual void setSeed( unsigned int seed );
	virtual unsigned int getSeed();
	virtual double getUniform();
	virtual TGenPhaseSpace* getGenerator() { return &event; }
};

class UniformTwoBodyGenerator : public RootGenerator
{
public:
	UniformTwoBodyGenerator(double minSq_, double maxSq_, int seed=-1) :
		RootGenerator(seed), minSq(minSq_),maxSq(maxSq_){
		if(Kinematics::instance()->getNumberOfParticles()!=2)
			throw std::runtime_error(
					"UniformTwoBodyGenerator::UniformTwoBodyGenerator() | Not a two body decay!");
	}
	virtual void generate(Event& evt);
	virtual UniformTwoBodyGenerator* Clone() { return (new UniformTwoBodyGenerator(*this)); }

protected:
	double minSq, maxSq;

};

#endif /* ROOTGENERATOR_HPP_ */

