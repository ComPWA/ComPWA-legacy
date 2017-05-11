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

namespace ComPWA {

class Event {

public:
  Event();

  Event(const std::string &name);

  Event(const double inWeight, const std::string &name,
        const double inEff = 1.);

  Event(Event const&);

  virtual void AddParticle(Particle inParticle);

  virtual ~Event();

  virtual void inline SetName(const std::string &name) { fName = name; }
  virtual const inline std::string &GetName() const { return fName; }
  virtual double inline GetWeight() const { return fWeight; }
  virtual void inline SetWeight(double w) { fWeight = w; }
  virtual int inline GetFlavour() const { return fFlavour; }
  virtual void inline SetFlavour(int fl) { fFlavour = fl; }
  virtual int inline GetCharge() const { return fCharge; }
  virtual void inline SetCharge(int ch) { fCharge = ch; }
  virtual double inline GetEfficiency() const { return fEff; }
  virtual void inline SetEfficiency(double eff) { fEff = eff; }

  virtual const inline unsigned long GetNParticles() const {
    return fParticles.size();
  }
  virtual const Particle& GetParticle(const unsigned int id) const;

  friend std::ostream &operator<<(std::ostream &stream, const Event &ev);

  virtual double GetCMSEnergy() const ;

  virtual void Clear() {
    fParticles.clear();
    fWeight = 1.0;
    fEff = 1.0;
    fName = "";
  }
protected:
  std::vector<Particle> fParticles;
  double fWeight;
  double fEff;
  std::string fName;
  int fFlavour; // 1 -> particle, 0 -> unknown, -1 anti-particle
  int fCharge;
};

} /* namespace ComPWA */

#endif
