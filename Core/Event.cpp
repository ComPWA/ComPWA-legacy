// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <algorithm>

#include "Core/Event.hpp"

namespace ComPWA {

Event::Event() : Weight(1.0) {}

std::ostream &operator<<(std::ostream &os, const Event &ev) {
  os << "Event: weight=" << ev.Weight << std::endl;
  os << " Printing particles (N=" << ev.ParticleList.size()
     << "):" << std::endl;
  for (auto const &x : ev.ParticleList)
    os << x << std::endl;

  return os;
}

double calculateInvariantMass(const Event &ev) {
  FourMomentum p4;
  for (auto x : ev.ParticleList)
    p4 += x.fourMomentum();
  return p4.invMass();
}

double getMaximumSampleWeight(const std::vector<Event> &sample) {
  double MaxWeight(0.0);
  auto MaxIterator = std::max_element(
      sample.begin(), sample.end(), [](const Event &a, const Event &b) -> bool {
        return a.Weight < b.Weight;
      });
  if (MaxIterator != sample.end())
    MaxWeight = MaxIterator->Weight;
  return MaxWeight;
}

DataPoint::DataPoint() : Weight(1.0) {}

std::ostream &operator<<(std::ostream &os, const DataPoint &p) {
  os << "(";
  for (unsigned int i = 0; i < p.KinematicVariableList.size(); ++i) {
    os << p.KinematicVariableList[i];
    if (i != p.KinematicVariableList.size() - 1) {
      os << ", ";
    }
  }
  os << ")";

  return os;
}

} // namespace ComPWA
