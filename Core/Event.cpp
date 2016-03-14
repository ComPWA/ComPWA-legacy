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

Event::Event():fWeight(1.),fName(""){

}

Event::Event(const std::string& name):fWeight(1.),fName(name),fFlavour(0),fCharge(0){

}

Event::Event(const double inWeight, const std::string& name="", const double inEff):fWeight(inWeight),fName(name),fFlavour(0),fCharge(0),fEff(inEff){

}

void Event::addParticle(Particle inParticle){
  fParticles.push_back(inParticle);
}

Event::~Event() { /* nothing */	}

const Particle& Event::getParticle(const unsigned int id) const{
  if(id>=getNParticles()){
    //TODO Exception
    return Particle();
  }
  return fParticles.at(id);
}

} /* namespace ComPWA */
