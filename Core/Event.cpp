// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <vector>
#include <string>
#include "Core/Particle.hpp"
#include "Core/Event.hpp"

using namespace ComPWA;

Event::Event() : Weight(1.), Eff(1.), Charge(0) {}

void Event::addParticle(Particle inParticle) {
  Particles.push_back(inParticle);
}

Particle Event::particle(size_t id) const {
  if (id >= numParticles()) {
    throw std::runtime_error(
        "Event::getParticle() | Particle id does not match!");
  }
  return Particles.at(id);
}

std::ostream &ComPWA::operator<<(std::ostream &stream, const Event &ev) {
  stream << "Event: weight=" << ev.weight() << " efficiency=" << ev.efficiency()
         << " charge=" << ev.charge() << std::endl;
  stream << " Printing particles (N=" << ev.numParticles() << "):" << std::endl;
  for (unsigned int i = 0; i < ev.numParticles(); ++i)
    stream << ev.particle(i) << std::endl;

  return stream;
}

double Event::cmsEnergy() const {

  FourMomentum p4;
  for (auto i : Particles)
    p4 += i.fourMomentum();
  return p4.invMass();
}

void Event::clear() {
  Particles.clear();
  Weight = 1.0;
  Eff = 1.0;
}
