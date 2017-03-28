#include "Core/TwoBodyKinematics.hpp"
#include "Core/PhysConst.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

TwoBodyKinematics::TwoBodyKinematics(int idMother, std::vector<int> finalState,
                                     double deltaMassWindow)
  : Kinematics(std::vector<int>(idMother), finalState) {
    
    assert( finalState.size() == 2 );
  _M = ComPWA::PhysConst::Instance()->FindParticle(idMother).GetMass();
  m1 = ComPWA::PhysConst::Instance()->FindParticle(finalState.at(0)).GetMass();
  m2 = ComPWA::PhysConst::Instance()->FindParticle(finalState.at(1)).GetMass();

  _spinM = ComPWA::PhysConst::Instance()
      ->FindParticle(idMother).GetSpin();
  spin1 = ComPWA::PhysConst::Instance()
      ->FindParticle(finalState.at(0)).GetSpin();
    spin2 = ComPWA::PhysConst::Instance()
      ->FindParticle(finalState.at(1)).GetSpin();
 

  if (_M == -999 || m1 == -999 || m2 == -999)
    throw std::runtime_error("TwoBodyKinematics(): Masses not set!");

  mass_min = ((_M - deltaMassWindow));
  mass_max = ((_M + deltaMassWindow));
  mass_sq_max = mass_max * mass_max;
  mass_sq_min = mass_min * mass_max;

  init();
}

void TwoBodyKinematics::init() {}

bool TwoBodyKinematics::IsWithinPhsp(const dataPoint &point) const {
  return 1;
  if (point.GetValue(0) >= mass_sq_min && point.GetValue(0) <= mass_sq_max)
    return 1;
  return 0;
}

void TwoBodyKinematics::EventToDataPoint(const Event &ev,
                                         dataPoint &point) const {
  double weight = ev.getWeight();
  point.SetWeight(weight); // reset weight
  const Particle &part1 = ev.getParticle(0);
  const Particle &part2 = ev.getParticle(1);
  double msq = Particle::invariantMass(part1, part2);
  point.SetValue(0, msq);
  return;
}

//! get mass of particles
double TwoBodyKinematics::GetMass(unsigned int num) const {
  if (num == 0)
    return _M;
  if (num == 1)
    return m1;
  if (num == 2)
    return m2;
  throw std::runtime_error("TwoBodyKinematics::getMass(int) | "
                           "Wrong particle requested!");
  return -999;
}

} /* namespace ComPWA */
