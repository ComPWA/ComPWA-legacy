// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <vector>
#include <string>
#include "Core/Particle.hpp"
#include "Core/Event.hpp"

namespace ComPWA {
Event::Event() : fWeight(1.), fEff(1.), fName(""), fFlavour(0), fCharge(0) {}

Event::Event(const std::string &name)
    : fWeight(1.), fEff(1.), fName(name), fFlavour(0), fCharge(0) {}

Event::Event(const double inWeight, const std::string &name = "",
             const double inEff)
    : fWeight(inWeight), fEff(inEff), fName(name), fFlavour(0), fCharge(0) {}

Event::Event(const Event &in)
    : fParticles(in.fParticles), fWeight(in.GetWeight()),
      fEff(in.GetEfficiency()), fName(in.GetName()), fFlavour(in.GetFlavour()),
      fCharge(in.GetCharge()) {}

void Event::AddParticle(Particle inParticle) {
  fParticles.push_back(inParticle);
}

Event::~Event() { /* nothing */
}

const Particle &Event::GetParticle(const unsigned int id) const {
  if (id >= GetNParticles()) {
    throw std::runtime_error(
        "Event::getParticle() | Particle id does not match!");
  }
  return fParticles.at(id);
}

std::ostream &operator<<(std::ostream &stream, const Event &ev) {
  stream << "Event name=" << ev.fName << " weight=" << ev.fWeight
         << " efficiency=" << ev.fEff << " flavour=" << ev.fFlavour
         << " charge=" << ev.fCharge << std::endl;
  stream << " Printing particles (N=" << ev.fParticles.size()
         << "):" << std::endl;
  for (auto i : ev.fParticles)
    stream << i << std::endl;

  return stream;
}

double Event::GetCMSEnergy() const {

  FourMomentum p4;
  for (auto i : fParticles)
    p4 += i.GetFourMomentum();
  return p4.GetInvMass();
}

} /* namespace ComPWA */
