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
private:
	TGenPhaseSpace event;
public:
	RootGenerator(int seed=-1){
		gRandom = new TRandom3(0);
		if(seed!=-1) setSeed(seed);
		DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
		//Generation
		TLorentzVector W(0.0, 0.0, 0.0, kin->getMass(0));//= beam + target;
		//(Momentum, Energy units are Gev/C, GeV)
		Double_t masses[3] = { kin->getMass(1), kin->getMass(2), kin->getMass(3)} ;
		event.SetDecay(W, 3, masses);
	};
	~RootGenerator(){};

	RootGenerator(const RootGenerator& other):event(other.event){}

	virtual RootGenerator* Clone() {
		return (new RootGenerator(*this));
	}

	virtual void generate(Event& evt) {
		TLorentzVector *p1,*p2,*p3;
		const double weight = event.Generate();

		p1 = event.GetDecay(0);
		p2 = event.GetDecay(1);
		p3 = event.GetDecay(2);

		Event tmp(weight,"");
		tmp.addParticle(Particle(p1->X(), p1->Y(), p1->Z(), p1->E()));
		tmp.addParticle(Particle(p2->X(), p2->Y(), p2->Z(), p2->E()));
		tmp.addParticle(Particle(p3->X(), p3->Y(), p3->Z(), p3->E()));
		evt=tmp;
	}

	void setSeed( unsigned int seed ){
		gRandom->SetSeed(seed);
	}

	unsigned int getSeed() { return gRandom->GetSeed(); }
	double getUniform(){
		return gRandom->Uniform(0,1);
	}
};


#endif /* ROOTGENERATOR_HPP_ */

