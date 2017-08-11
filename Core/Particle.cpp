// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Particle.hpp"

namespace ComPWA {

Particle::Particle(double inPx, double inPy, double inPz, double inE, int inpid,
                   int c)
  :  _p4(std::array<double,4>{{inPx, inPy, inPz, inE}}), pid(inpid), charge(c)
  {
  }

Particle::Particle(const Particle &in)
    : _p4(in._p4), pid(in.pid), charge(in.charge) {}

Particle::~Particle() {}

std::ostream &operator<<(std::ostream &stream, const Particle &p) {
  stream << "Particle id=" << p.pid << " charge=" << p.charge
         << " p4=" << p.GetFourMomentum();
  return stream;
}

// Particle &Particle::operator+=(const Particle &rhs) {
//  std::transform(_p4.begin(), _p4.end(), rhs._p4.begin(), _p4.begin(),
//                 std::plus<double>());
//  charge += rhs.GetCharge();
//  pid += rhs.GetPid(); // doesn't really make sense here
//  return *this;
//}
//
// Particle Particle::operator+(const Particle &c2) {
//  Particle newP(*this);
//  newP += c2;
//  return newP;
//}

double Particle::InvariantMass(const Particle &inPa, const Particle &inPb) {
  return FourMomentum::InvariantMass(inPa.GetFourMomentum(), inPb.GetFourMomentum());
}

} /* namespace ComPWA */
