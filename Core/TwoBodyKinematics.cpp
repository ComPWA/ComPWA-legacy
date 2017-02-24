#include "Core/TwoBodyKinematics.hpp"
#include "Core/PhysConst.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

TwoBodyKinematics::TwoBodyKinematics(std::string _nameMother,
                                     std::string _name1, std::string _name2,
                                     double deltaMassWindow)
    : Kinematics(_nameMother, 0.0, 2), name1(_name1), name2(_name2) {
  _M = ComPWA::PhysConst::instance()->findParticle(_nameMother).GetMass();
  m1 = ComPWA::PhysConst::instance()->findParticle(_name1).GetMass();
  m2 = ComPWA::PhysConst::instance()->findParticle(_name2).GetMass();

  _spinM = ComPWA::PhysConst::instance()
      ->findParticle(_nameMother).GetSpin();
  spin1 = ComPWA::PhysConst::instance()
      ->findParticle(_name1).GetSpin();
    spin2 = ComPWA::PhysConst::instance()
      ->findParticle(_name2).GetSpin();
 

  if (_M == -999 || m1 == -999 || m2 == -999)
    throw std::runtime_error("TwoBodyKinematics(): Masses not set!");

  mass_min = ((_M - deltaMassWindow));
  mass_max = ((_M + deltaMassWindow));
  mass_sq_max = mass_max * mass_max;
  mass_sq_min = mass_min * mass_max;
  _varNames.push_back("msq");

  init();
}

void TwoBodyKinematics::init() {}

bool TwoBodyKinematics::IsWithinPhsp(const dataPoint &point) const {
  return 1;
  if (point.getVal(0) >= mass_sq_min && point.getVal(0) <= mass_sq_max)
    return 1;
  return 0;
}

void TwoBodyKinematics::EventToDataPoint(const Event &ev,
                                         dataPoint &point) const {
  double weight = ev.getWeight();
  point.setWeight(weight); // reset weight
  const Particle &part1 = ev.getParticle(0);
  const Particle &part2 = ev.getParticle(1);
  double msq = Particle::invariantMass(part1, part2);
  point.setVal(0, msq);
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

//! get mass of paticles
double TwoBodyKinematics::GetMass(std::string name) const {
  if (name == _nameMother)
    return _M;
  if (name == name1)
    return m1;
  if (name == name2)
    return m2;
  throw std::runtime_error("TwoBodyKinematics::getMass(int) | "
                           "Wrong particle " +
                           name + " requested!");
  return -999;
}

} /* namespace ComPWA */
