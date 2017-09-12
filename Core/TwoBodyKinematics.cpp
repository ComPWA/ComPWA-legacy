// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/TwoBodyKinematics.hpp"
#include "Core/Properties.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

TwoBodyKinematics::TwoBodyKinematics(std::shared_ptr<PartList> partL,
                                     int idMother, std::vector<int> finalState,
                                     double deltaMassWindow)
    : Kinematics(std::vector<int>(idMother), finalState) {

  assert(finalState.size() == 2);
  auto propM = FindParticle(partL, idMother);
  auto propM1 = FindParticle(partL, finalState.at(0));
  auto propM2 = FindParticle(partL, finalState.at(1));
  _M = propM.GetMass();
  m1 = propM1.GetMass();
  m2 = propM2.GetMass();

  _spinM = propM.GetSpinQuantumNumber("Spin");

  if (_M == -999 || m1 == -999 || m2 == -999)
    throw std::runtime_error("TwoBodyKinematics(): Masses not set!");

  mass_min = ((_M - deltaMassWindow));
  mass_max = ((_M + deltaMassWindow));
  mass_sq_max = mass_max * mass_max;
  mass_sq_min = mass_min * mass_max;

}

bool TwoBodyKinematics::IsWithinPhsp(const dataPoint &point) const {
  if (point.GetValue(0) >= mass_sq_min && point.GetValue(0) <= mass_sq_max)
    return 1;
  return 0;
}

void TwoBodyKinematics::EventToDataPoint(const Event &ev,
                                         dataPoint &point) const {
  double weight = ev.GetWeight();
  point.SetWeight(weight); // reset weight
  const Particle &part1 = ev.GetParticle(0);
  const Particle &part2 = ev.GetParticle(1);
  double msq = Particle::InvariantMass(part1, part2);
  point.SetValue(0, msq);
  return;
}

//! get mass of particles
double TwoBodyKinematics::GetMass(unsigned int num) const {
  if( num > 2 )
  throw std::runtime_error("TwoBodyKinematics::getMass(int) | "
                           "Wrong particle requested!");
  if (num == 0)
    return _M;
  if (num == 1)
    return m1;
  
  return m2;
}

} /* namespace ComPWA */
