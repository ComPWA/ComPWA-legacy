// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Event class
///

#ifndef _Event_HPP_
#define _Event_HPP_

#include <vector>
#include <string>
#include <iostream>
#include "Core/Particle.hpp"

namespace ComPWA {

///
/// \class Event
/// Internal container for event information. This class provides a internal
/// container for event-based information. The class provides a list of
/// particles of the event.
///
class Event {

public:
  Event();

  Event(const std::string &name);

  Event(const double inWeight, const std::string &name,
        const double inEff = 1.);

  virtual void addParticle(Particle inParticle);

  virtual ~Event(){};

  virtual double inline weight() const { return Weight; }

  virtual void inline setWeight(double w) { Weight = w; }

  virtual int inline charge() const { return Charge; }

  virtual void inline setCharge(int ch) { Charge = ch; }

  virtual double inline efficiency() const { return Eff; }

  virtual void inline setEfficiency(double eff) { Eff = eff; }

  virtual size_t numParticles() const { return Particles.size(); }

  virtual Particle particle(size_t id) const;

  virtual std::vector<Particle> &particles() { return Particles; }

  friend std::ostream &operator<<(std::ostream &stream, const Event &ev);

  /// Invariant mass of all particles in the event
  virtual double cmsEnergy() const;

  virtual void clear();

  std::vector<Particle>::iterator first() { return Particles.begin(); }

  std::vector<Particle>::iterator last() { return Particles.end(); }

protected:
  std::vector<Particle> Particles;
  double Weight;
  double Eff;
  int Charge;
};

} // namespace ComPWA

#endif
