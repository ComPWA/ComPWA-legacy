// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <algorithm>

#include "Core/Event.hpp"

namespace ComPWA {

std::ostream &operator<<(std::ostream &os, const Event &ev) {
  os << "Event: weight=" << ev.Weight << "\n";
  os << "Particle four-momenta (N=" << ev.FourMomenta.size() << "):\n";
  for (auto const &x : ev.FourMomenta)
    os << x << std::endl;

  return os;
}

double calculateInvariantMass(const Event &ev) {
  FourMomentum p4;
  for (auto x : ev.FourMomenta)
    p4 += x;
  return p4.invariantMass();
}

double getMaximumSampleWeight(const EventCollection &Sample) {
  double MaxWeight(0.0);
  auto MaxIterator =
      std::max_element(Sample.Events.begin(), Sample.Events.end(),
                       [](const Event &a, const Event &b) -> bool {
                         return a.Weight < b.Weight;
                       });
  if (MaxIterator != Sample.Events.end())
    MaxWeight = MaxIterator->Weight;
  return MaxWeight;
}

} // namespace ComPWA
