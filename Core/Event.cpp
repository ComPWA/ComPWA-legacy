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
//-------------------------------------------------------------------------------
#include <vector>
#include <string>
#include "Core/Particle.hpp"
#include "Core/Event.hpp"

namespace ComPWA {
Event::Event() : fWeight(1.),fEff(1.),fName(""),fFlavour(0),fCharge(0)
{

}

Event::Event(const std::string& name):fWeight(1.),fEff(1.),fName(name),fFlavour(0),fCharge(0){

}

Event::Event(const double inWeight, const std::string& name="", const double inEff):fWeight(inWeight),fEff(inEff),fName(name),fFlavour(0),fCharge(0){

}

void Event::addParticle(Particle inParticle){
  fParticles.push_back(inParticle);
}

void Event::setParticleAt(const Particle &particle, unsigned int index) {
  fParticles[index] = particle;
}

Event::~Event() { /* nothing */	}

const Particle& Event::getParticle(const unsigned int id) const{
  if(id>=getNParticles()){
      throw std::runtime_error("Event::getParticle() | Particle id does not match!");
  }
  return fParticles.at(id);
}

std::ostream& operator<< (std::ostream& stream, const Event& ev){
	stream<< "Event name="<<ev.fName<<" weight="<<ev.fWeight<<" efficiency="<<ev.fEff
			<<" flavour="<<ev.fFlavour<<" charge="<<ev.fCharge<<std::endl;
	stream<< " Printing particles (N="<<ev.fParticles.size()<<"):"<<std::endl;
	for(auto i: ev.fParticles)
		stream<<i<<std::endl;

	return stream;
}
} /* namespace ComPWA */
