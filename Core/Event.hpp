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
//! Internal container for event information.
/*! \class Event
 * @file Event.hpp
 * This class provides a internal container for event-based information. The
 * class provides a list of particles of the event.
*/

#ifndef _Event_HPP_
#define _Event_HPP_

#include <vector>
#include <string>
#include <iostream>
#include "Core/Particle.hpp"

class Event{

public:
  Event();

  Event(const std::string& name);

  Event(const double inWeight, const std::string& name, const double inEff = 1.);

  virtual void addParticle(Particle inParticle);

  virtual ~Event();

  virtual void inline setName(const std::string& name) { fName = name; }
  virtual const inline std::string& getName() const { return fName; }
  virtual double inline getWeight() const {return fWeight;}
  virtual void inline setWeight(double w) { fWeight=w;}
  virtual int inline getFlavour() const {return fFlavour;}
  virtual void inline setFlavour(int fl) { fFlavour = fl;}
  virtual int inline getCharge() const {return fCharge;}
  virtual void inline setCharge(int ch) { fCharge = ch;}
  virtual double inline getEfficiency() const {return fEff;}
  virtual void inline setEfficiency(double eff) { fEff = eff;}

  virtual const inline unsigned int getNParticles() const { return fParticles.size(); }
  virtual const Particle& getParticle(const unsigned int id);

  friend std::ostream& operator<< (std::ostream& stream, const Event& ev);

protected:
  std::vector<Particle> fParticles;
  double fWeight;
  double fEff;
  std::string fName;
  int fFlavour; //1 -> particle, 0 -> unknown, -1 anti-particle
  int fCharge;

};

#endif
