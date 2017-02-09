/*
 * Particle.cpp
 *
 *  Created on: Nov 16, 2015
 *      Author: weidenka
 */

#include "Core/Particle.hpp"

namespace ComPWA {

Particle::Particle(double inPx, double inPy, double inPz, double inE, int inpid,
                   int c)
    : px(inPx), py(inPy), pz(inPz), E(inE), pid(inpid), charge(c) {}

Particle::Particle(const Particle &in)
    : px(in.getPx()), py(in.getPy()), pz(in.getPz()), E(in.getE())
    , pid(in.getPid()), charge(in.getCharge()) {}

Particle::~Particle() {}

std::ostream &operator<<(std::ostream &stream, const Particle &p) {
  stream << "Particle id=" << p.pid << " charge=" << p.charge << " p4=(" << p.px
         << "," << p.py << "," << p.pz << "," << p.px << ")";
  return stream;
}

Particle &Particle::operator+=(const Particle &rhs) {
  px += rhs.getPx();
  py += rhs.getPy();
  pz += rhs.getPz();
  E += rhs.getE();
  charge += rhs.getCharge();
  pid += rhs.getPid(); // doesn't really make sense here
  return *this;
}

Particle Particle::operator+(const Particle &c2) {
  // use the Cents constructor and operator+(int, int)
  Particle newP(*this);
  newP += c2;
  return newP;
}

double Particle::invariantMass(const Particle &inPa, const Particle &inPb) {
  return (std::pow(inPa.getE() + inPb.getE(), 2) -
          std::pow(inPa.getPx() + inPb.getPx(), 2) -
          std::pow(inPa.getPy() + inPb.getPy(), 2) -
          std::pow(inPa.getPz() + inPb.getPz(), 2));
}
}
