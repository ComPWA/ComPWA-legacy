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
  FourMomentum(double px = 0, double py = 0, double pz = 0, double E = 0)
      : P4(std::array<double, 4>{{px, py, pz, E}}) {}

  FourMomentum(std::array<double, 4> p4) : P4(p4) {}

  FourMomentum(std::vector<double> p4) {
    if (p4.size() != 4)
      throw std::runtime_error(
          "FourMomentum::Fourmomentum() | Size of vector not equal 4!");
    P4 = std::array<double, 4>{{p4.at(0), p4.at(1), p4.at(2), p4.at(3)}};
  }

  /*void setPx(double px) { P4.at(0) = px; }
  void setPy(double py) { P4.at(1) = py; }
  void setPz(double pz) { P4.at(2) = pz; }
  void setE(double E) { P4.at(3) = E; }*/
  double px() const { return P4.at(0); }
  double py() const { return P4.at(1); }
  double pz() const { return P4.at(2); }
  double e() const { return P4.at(3); }

  FourMomentum operator+(const FourMomentum &pB) const {
    FourMomentum newP(*this);
    newP += pB;
    return newP;
  }

  void operator+=(const FourMomentum &pB) {
    std::transform(P4.begin(), P4.end(), pB.P4.begin(), P4.begin(),
                   std::plus<double>());
  }

  operator std::vector<double>() {
    return std::vector<double>(P4.begin(), P4.end());
  }

  operator std::array<double, 4>() { return P4; }

  std::array<double, 4> operator()() const { return P4; }

  bool operator==(const FourMomentum &pB) const {
    if (P4 == pB.P4)
      return true;
    return false;
  }

  friend std::ostream &operator<<(std::ostream &stream,
                                  const FourMomentum &p4) {
    auto vec = p4.value();
    stream << "(" << vec.at(0) << "," << vec.at(1) << "," << vec.at(2) << ","
           << vec.at(3) << ")";
    return stream;
  }

  const std::array<double, 4> &value() const { return P4; }

  // void setValue(std::array<double, 4> p4) { P4 = p4; }

  double invMassSq() const { return invariantMass(*this); }

  double invMass() const { return std::sqrt(invMassSq()); }

  static double invariantMass(const FourMomentum &p4A,
                              const FourMomentum &p4B) {
    return invariantMass(p4A + p4B);
  }

  static double invariantMass(const FourMomentum &p4) {
    auto vec = p4.value();
    return ((-1) * (vec.at(0) * vec.at(0) + vec.at(1) * vec.at(1) +
                    vec.at(2) * vec.at(2) - vec.at(3) * vec.at(3)));
  }

  static double threeMomentumSq(const FourMomentum &p4) {
    auto vec = p4.value();
    return (vec.at(0) * vec.at(0) + vec.at(1) * vec.at(1) +
            vec.at(2) * vec.at(2));
  }

private:
  std::array<double, 4> P4;
};

///
/// \class Particle
/// ComPWA particle class.
/// This class provides a internal container for information of a particle. The
/// class provides the momentum 4-vector and pid of the particle.
///
class Particle {
public:
  Particle(double inPx = 0, double inPy = 0, double inPz = 0, double inE = 0,
           int inpid = 0);

  Particle(std::array<double, 4> p4, int inpid = 0) : P4(p4), Pid(inpid){};

  Particle(Particle const &);

  virtual ~Particle(){};

  /*virtual void px(double px) { P4.setPx(px); }
  virtual void py(double py) { P4.setPy(py); }
  virtual void pz(double pz) { P4.setPz(pz); }
  virtual void e(double E) { P4.setE(E); }
  virtual double px() const { return P4.px(); }
  virtual double py() const { return P4.py(); }
  virtual double pz() const { return P4.pz(); }
  virtual double e() const { return P4.e(); }*/

  // virtual void setPid(int _pid) { Pid = _pid; }
  virtual int pid() const { return Pid; }

  virtual const FourMomentum &fourMomentum() const { return P4; }

  friend std::ostream &operator<<(std::ostream &stream, const Particle &p) {
    stream << "Particle id=" << p.pid() << " p4=" << p.fourMomentum();
    return stream;
  }

  /// Magnitude of three momentum
  double threeMomentum() const {
    return std::sqrt(FourMomentum::threeMomentumSq(P4));
  }

  /// Get invariant mass
  virtual double mass() const { return std::sqrt(massSq()); }

  virtual double massSq() const { return FourMomentum::invariantMass(P4); }

  /// Invariant mass of \p inPa and \p inPb.
  static double invariantMass(const Particle &inPa, const Particle &inPb);

private:
  FourMomentum P4;
  int Pid;
};

} // namespace ComPWA
#endif
