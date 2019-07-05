// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Particle.hpp"

using namespace ComPWA;

Particle::Particle(double inPx, double inPy, double inPz, double inE, int inpid)
    : P4(std::array<double, 4>{{inPx, inPy, inPz, inE}}), Pid(inpid) {}

Particle::Particle(const Particle &in) : P4(in.P4), Pid(in.Pid) {}

double Particle::invariantMass(const Particle &inPa, const Particle &inPb) {
  return FourMomentum::invariantMass(inPa.fourMomentum(), inPb.fourMomentum());
}
