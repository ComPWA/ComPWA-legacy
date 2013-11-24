/*
 * RootGenerator.hpp
 *
 *  Created on: Nov 23, 2013
 *      Author: weidenka
 */

#ifndef ROOTGENERATOR_HPP_
#define ROOTGENERATOR_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include "Core/Generator.hpp"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"

class RootGenerator : public Generator {
private:
	TGenPhaseSpace event;
public:
	RootGenerator(){
		DalitzKinematics* kin = DalitzKinematics::instance();
		//Generation
		TLorentzVector W(0.0, 0.0, 0.0, kin->getMass(0));//= beam + target;
		//(Momentum, Energy units are Gev/C, GeV)
		Double_t masses[3] = { kin->getMass(1), kin->getMass(2), kin->getMass(3)} ;
		event.SetDecay(W, 3, masses);
	};
	~RootGenerator(){};
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
};


#endif /* ROOTGENERATOR_HPP_ */

