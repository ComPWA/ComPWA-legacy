// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Particle and FourMomentum class.
///

#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <boost/property_tree/ptree.hpp>

namespace ComPWA {

///
/// \class FourMomentum
/// ComPWA four momentum class.
///
class FourMomentum {

public:
  FourMomentum() : FourMomentum(0.0, 0.0, 0.0, 0.0){};
  FourMomentum(double E) : FourMomentum(0.0, 0.0, 0.0, E){};

  FourMomentum(double Px, double Py, double Pz, double E)
      : FourMomentum(std::array<double, 4>{{Px, Py, Pz, E}}) {}

  // main constructor
  FourMomentum(std::array<double, 4> P4_) : P4(P4_) {}

  FourMomentum(std::vector<double> P4_)
      : FourMomentum([&P4_]() {
          if (P4_.size() != 4)
            throw std::runtime_error(
                "FourMomentum::Fourmomentum() | Size of vector not equal 4!");
          return std::array<double, 4>{{P4_[0], P4_[1], P4_[2], P4_[3]}};
        }()) {}

  double px() const { return P4[0]; }
  double py() const { return P4[1]; }
  double pz() const { return P4[2]; }
  double e() const { return P4[3]; }

  FourMomentum operator+(const FourMomentum &pB) const {
    FourMomentum newP(*this);
    newP += pB;
    return newP;
  }

  void operator+=(const FourMomentum &pB) {
    std::transform(P4.begin(), P4.end(), pB.P4.begin(), P4.begin(),
                   std::plus<double>());
  }

  bool operator==(const FourMomentum &pB) const { return P4 == pB.P4; }

  friend std::ostream &operator<<(std::ostream &stream,
                                  const FourMomentum &p4) {
    stream << "(" << p4.px() << "," << p4.py() << "," << p4.py() << ","
           << p4.e() << ")";
    return stream;
  }

  double invariantMassSquared() const {
    return ((-1) *
            (P4[0] * P4[0] + P4[1] * P4[1] + P4[2] * P4[2] - P4[3] * P4[3]));
  }

  double invariantMass() const { return std::sqrt(invariantMassSquared()); }

  double threeMomentumSquared() const {
    return (P4[0] * P4[0] + P4[1] * P4[1] + P4[2] * P4[2]);
  }

private:
  std::array<double, 4> P4;
}; // namespace ComPWA

///
/// \class Particle
/// ComPWA particle class.
/// This class provides a internal container for information of a particle. The
/// class provides the momentum 4-vector and pid of the particle.
///
class Particle {
public:
  Particle(double inPx, double inPy, double inPz, double inE, int inpid)
      : Particle(std::array<double, 4>{{inPx, inPy, inPz, inE}}, inpid) {}

  Particle(std::array<double, 4> p4, int inpid) : P4(p4), Pid(inpid){};

  virtual ~Particle(){};

  virtual int pid() const { return Pid; }

  virtual const FourMomentum &fourMomentum() const { return P4; }

  friend std::ostream &operator<<(std::ostream &stream, const Particle &p) {
    stream << "Particle id=" << p.pid() << " p4=" << p.fourMomentum();
    return stream;
  }

  /// Magnitude of three momentum
  double threeMomentum() const { return std::sqrt(P4.threeMomentumSquared()); }

  /// Get invariant mass
  virtual double mass() const { return P4.invariantMass(); }

  virtual double massSquared() const { return P4.invariantMassSquared(); }

  /// Invariant mass of \p inPa and \p inPb.
  static double invariantMass(const Particle &inPa, const Particle &inPb) {
    return (inPa.fourMomentum() + inPb.fourMomentum()).invariantMass();
  }

private:
  FourMomentum P4;
  int Pid;
};

} // namespace ComPWA
#endif
