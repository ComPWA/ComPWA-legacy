//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Internal container representing a particle.
/*! \class Particle
 * @file Particle.hpp
 * This class provides a internal container for information of a particle. The
 * class provides the momentum 4-vector and pid of the particle.
 */

#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include <vector>
#include <array>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace ComPWA {

class FourMomentum {

public:
  FourMomentum(double px = 0, double py = 0, double pz = 0, double E = 0)
      : _p4(std::array<double, 4>{{px, py, pz, E}}) {
    //    _p4.at(0) = px;
    //    _p4.at(1) = py;
    //    _p4.at(2) = pz;
    //    _p4.at(3) = E;
  }
  FourMomentum(std::array<double, 4> p4) : _p4(p4) {}

  void SetPx(double px) { _p4.at(0) = px; }
  void SetPy(double py) { _p4.at(1) = py; }
  void SetPz(double pz) { _p4.at(2) = pz; }
  void SetE(double E) { _p4.at(3) = E; }
  double GetPx() const { return _p4.at(0); }
  double GetPy() const { return _p4.at(1); }
  double GetPz() const { return _p4.at(2); }
  double GetE() const { return _p4.at(3); }

  FourMomentum operator+(const FourMomentum &pB) const {
    FourMomentum newP(*this);
    newP += pB;
    return newP;
  }

  void operator+=(const FourMomentum &pB) {
    std::transform(_p4.begin(), _p4.end(), pB._p4.begin(), _p4.begin(),
                   std::plus<double>());
  }
  bool operator==(const FourMomentum &pB) const {
    if (_p4 == pB._p4)
      return true;
    return false;
  }

  friend std::ostream &operator<<(std::ostream &stream,
                                  const FourMomentum &p4) {
    stream << "(" << p4.GetPx() << "," << p4.GetPy() << "," << p4.GetPz() << ","
           << p4.GetE() << ")";
    return stream;
  }

  virtual inline const std::array<double, 4> &GetFourMomentum() const {
    return _p4;
  }

  virtual void SetFourMomentum(std::array<double, 4> p4) { _p4 = p4; }

  inline double GetMassSq() const { return InvariantMass(*this); }

  static double InvariantMass(const FourMomentum &p4A,
                              const FourMomentum &p4B) {
    return InvariantMass(p4A + p4B);
  }

  static double InvariantMass(const FourMomentum &p4) {
    auto vec = p4.GetFourMomentum();
    return ((-1) * (vec.at(0) * vec.at(0) + vec.at(1) * vec.at(1) +
                    vec.at(2) * vec.at(2) - vec.at(3) * vec.at(3)));
  }

  static double ThreeMomentumSq(const FourMomentum &p4) {
    auto vec = p4.GetFourMomentum();
    return (vec.at(0) * vec.at(0) + vec.at(1) * vec.at(1) +
            vec.at(2) * vec.at(2));
  }

protected:
  std::array<double, 4> _p4;
};

class Particle {
public:
  //! Default constructor
  Particle(double inPx = 0, double inPy = 0, double inPz = 0, double inE = 0,
           int inpid = 0, int c = 0);

  Particle(std::array<double, 4> p4, int inpid = 0, int c = 0)
      : _p4(p4), pid(inpid), charge(c){};

  //! Copy constructor
  Particle(Particle const &);

  //! Default destructor
  virtual ~Particle();

  virtual inline void SetPx(double px) { _p4.SetPx(px); }
  virtual inline void SetPy(double py) { _p4.SetPy(py); }
  virtual inline void SetPz(double pz) { _p4.SetPz(pz); }
  virtual inline void SetE(double E) { _p4.SetE(E); }
  virtual inline double GetPx() const { return _p4.GetPx(); }
  virtual inline double GetPy() const { return _p4.GetPy(); }
  virtual inline double GetPz() const { return _p4.GetPz(); }
  virtual inline double GetE() const { return _p4.GetE(); }

  virtual inline void SetPid(int _pid) { pid = _pid; }
  virtual inline void SetCharge(int _c) { charge = _c; }
  virtual inline int GetPid() const { return pid; }
  virtual inline int GetCharge() const { return charge; }

  virtual inline const FourMomentum &GetFourMomentum() const { return _p4; }

  friend std::ostream &operator<<(std::ostream &stream, const Particle &p);

  //  Particle &operator+=(const Particle &rhs);

  //  Particle operator+(const Particle &c2);

  //! Get magnitude of three momentum
  double GetThreeMomentum() const {
    return std::sqrt(FourMomentum::ThreeMomentumSq(_p4));
  }

  //! Get invariant mass
  double GetMass() const { return std::sqrt(GetMassSq()); }

  inline double GetMassSq() const { return FourMomentum::InvariantMass(_p4); }

  /**! Invariant mass
   *
   * @param inPa Particle A
   * @param inPb Particle B
   * @return
   */
  static double InvariantMass(const Particle &inPa, const Particle &inPb);

protected:
  FourMomentum _p4;
  int pid;
  int charge;
};

} /* namespace ComPWA */
#endif
