// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Spin class.
///

#ifndef CORE_UTILITY_HPP_
#define CORE_UTILITY_HPP_

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include "Core/Exceptions.hpp"

namespace ComPWA {

typedef std::vector<unsigned int> IndexList;
typedef std::pair<unsigned int, unsigned int> IndexPair;
typedef std::map<unsigned int, unsigned int> IndexMapping;

///
/// \class Spin
/// Spin class for integer and half-integer spins.
///
class Spin {

public:
  Spin()
      : J_numerator_(0), J_denominator_(1), J_z_numerator_(0),
        z_component_relevant(false) {}

  /// Constructor for half-integer spin
  Spin(int num, int denom)
      : J_numerator_(num), J_denominator_(denom), J_z_numerator_(0),
        z_component_relevant(false) {}

  /// Constructor for integer spin
  Spin(unsigned int intSpin)
      : J_numerator_(intSpin), J_denominator_(1), J_z_numerator_(0),
        z_component_relevant(false) {}

  /// Constructor for integer spin
  Spin(int intSpin)
      : J_numerator_(intSpin), J_denominator_(1), J_z_numerator_(0),
        z_component_relevant(false) {}

  /// Constructor for double spin
  Spin(double spin, double spinZ = 0.0)
      : J_numerator_(0), J_denominator_(1), J_z_numerator_(0),
        z_component_relevant(false) {
        // Make sure all variables are initialized
    if (isInteger(spin)) {
      SetNumerator(spin);
      SetDenominator(1);
      if (isInteger(spinZ))
        SetZNumerator(spinZ);
    } else if (isInteger(2 * spin)) {
      SetNumerator(2 * spin);
      SetDenominator(2);
      if (isInteger(2 * spinZ))
        SetZNumerator(2 * spinZ);
    } else
      throw BadParameter("Spin::Spin() |" + std::to_string(spin) +
                         " is not a valid spin!");
  }

  bool equalMagnitude(const Spin &rhs) const {
    if (1.0 * this->J_numerator_ / this->J_denominator_ ==
        1.0 * rhs.J_numerator_ / rhs.J_denominator_)
      return true;
    return false;
  }

  Spin operator+(const Spin &rhs) const {
    // We calculate (a/b - c/d) = (a*d - c*b)/(bd)
    int num = GetNumerator() * rhs.GetDenominator() +
              rhs.GetNumerator() * GetDenominator();
    int denom = rhs.GetDenominator() * GetDenominator();
    int znum = GetZNumerator() * rhs.GetDenominator() +
               rhs.GetZNumerator() * GetDenominator();
    Spin s;
    s.SetNumerator(num);
    s.SetZNumerator(znum);
    s.SetDenominator(denom);
    s.Simplify();
    return s;
  }

  Spin operator-(const Spin &rhs) const {
    // We calculate (a/b - c/d) = (a*d - c*b)/(bd)
    int num = GetNumerator() * rhs.GetDenominator() -
              rhs.GetNumerator() * GetDenominator();
    int denom = rhs.GetDenominator() * GetDenominator();
    int znum = GetZNumerator() * rhs.GetDenominator() -
               rhs.GetZNumerator() * GetDenominator();
    Spin s;
    s.SetNumerator(num);
    s.SetZNumerator(znum);
    s.SetDenominator(denom);
    s.Simplify();
    return s;
  }

  Spin &operator=(const int &other) {
    SetNumerator(other);
    SetDenominator(1);
    return *this;
  }

  Spin &operator=(const double &other) {
    if (isInteger(other)) {
      SetNumerator(other);
      SetDenominator(1);
    } else if (isInteger(2 * other)) {
      SetNumerator(2 * other);
      SetDenominator(2);
    } else
      throw BadParameter("Spin::operator=() |" + std::to_string(other) +
                         " is not a valid spin!");

    return *this;
  }

  Spin &operator--() {
    // Only decrement if the result is larger zero
    if (((int)GetNumerator()) - 1 >= 0) {
      SetNumerator(GetNumerator() - 1);
      if (z_component_relevant)
        SetZNumerator(GetZNumerator() + 1);
    }
    return *this;
  }

  Spin operator--(int) {
    Spin tmp(*this); // copy
    operator++();    // pre-increment
    return tmp;      // return old value
  }

  Spin &operator++() {
    SetNumerator(GetNumerator() + 1);
    if (z_component_relevant)
      SetZNumerator(GetZNumerator() + 1);
    return *this;
  }

  Spin operator++(int) {
    Spin tmp(*this); // copy
    operator++();    // pre-increment
    return tmp;      // return old value
  }

  Spin operator*(const int factor) {
    SetNumerator(GetNumerator() * factor);
    return *this; // return new value
  }

  // conversion to double (type-cast operator)
  operator double() const { return ((double)J_numerator_) / J_denominator_; }

  void Simplify() {
    int tmp, tmpZ;
    tmp = tmpZ = ggT(GetNumerator(), GetDenominator());
    if (UseZ())
      tmpZ = ggT(GetZNumerator(), GetDenominator());
    if (tmp == tmpZ) {
      SetNumerator(GetNumerator() / tmp);
      SetZNumerator(GetZNumerator() / tmp);
      SetDenominator(GetDenominator() / tmp);
    }
  }

  /// Calculate largest common factor
  static unsigned int ggT(unsigned int a, unsigned int b) {
    if (b == 0)
      return a;
    else
      return ggT(b, a % b);
  }

  static unsigned int kgV(unsigned int a, unsigned int b) {
    return (a * b) / ggT(a, b);
  }

  static bool isInteger(double a) {
    if (a == (int)a)
      return true;
    return false;
  }

  bool operator==(const Spin &rhs) const {
    if (this->J_numerator_ != rhs.J_numerator_)
      return false;
    if (this->J_denominator_ != rhs.J_denominator_)
      return false;
    if (this->J_z_numerator_ != rhs.J_z_numerator_)
      return false;

    return true;
  }

  bool operator!=(const Spin &rhs) const { return !(*this == rhs); }

  bool operator<(const Spin &rhs) const {
    if (this->J_numerator_ < rhs.J_numerator_)
      return true;
    else if (this->J_numerator_ > rhs.J_numerator_)
      return false;
    if (this->J_denominator_ < rhs.J_denominator_)
      return true;
    else if (this->J_denominator_ > rhs.J_denominator_)
      return false;
    if (this->J_z_numerator_ < rhs.J_z_numerator_)
      return true;

    return false;
  }

  bool operator<=(const Spin &rhs) const {
    if (*this == rhs)
      return true;
    if (*this < rhs)
      return true;
    return false;
  }

  bool operator>(const Spin &rhs) const { return (rhs < *this); }

  bool operator>=(const Spin &rhs) {
    if (*this == rhs)
      return true;
    if (*this > rhs)
      return true;
    return false;
  }

  int GetNumerator() const { return J_numerator_; }

  int GetZNumerator() const { return J_z_numerator_; }

  unsigned int GetDenominator() const { return J_denominator_; }

  void SetNumerator(unsigned int num) { J_numerator_ = num; }

  void SetZNumerator(int znum) { J_z_numerator_ = znum; }

  void SetDenominator(unsigned int denom) {
    if (denom != 1 && denom != 2)
      throw BadParameter("Spin::SetDenominator() |"
                         " Should be equal 1 oder 2!");
    J_denominator_ = denom;
  }

  double GetSpin() { return (double)J_numerator_ / J_denominator_; }

  double GetZComponent() { return (double)J_z_numerator_ / J_denominator_; }

  bool UseZ() const { return z_component_relevant; }

  void SetUseZ(bool b) { z_component_relevant = b; }

protected:
  int J_numerator_;
  
  unsigned int J_denominator_;
  
  int J_z_numerator_;

  bool z_component_relevant;
};

struct IDInfo {
  int particleId_;
  std::string name_;

  bool operator==(const IDInfo &rhs) const {
    if (this->particleId_ != rhs.particleId_)
      return false;
    if (this->name_ != rhs.name_)
      return false;

    return true;
  }
  bool operator!=(const IDInfo &rhs) const { return !(*this == rhs); }

  bool operator<(const IDInfo &rhs) const {
    return lessThenIgnoringID(*this, rhs);
  }

  static bool lessThenIgnoringID(const IDInfo &lhs, const IDInfo &rhs) {
    if (lhs.particleId_ < rhs.particleId_)
      return true;
    else if (lhs.particleId_ > rhs.particleId_)
      return false;
    if (lhs.name_ < rhs.name_)
      return true;

    return false;
  }

  bool operator>(const IDInfo &rhs) const { return (rhs < *this); }
};

struct ParticleStateInfo {
  unsigned int uniqueId_;
  IDInfo pid_information_;
  Spin spin_information_;
  bool coherent;

  bool operator==(const ParticleStateInfo &rhs) const {
    if (this->uniqueId_ != rhs.uniqueId_)
      return false;
    if (this->pid_information_ != rhs.pid_information_)
      return false;
    if (this->spin_information_ != rhs.spin_information_)
      return false;
    if (this->coherent != rhs.coherent)
      return false;

    return true;
  }

  bool operator!=(const ParticleStateInfo &rhs) const {
    return !((*this) == rhs);
  }

  bool operator<(const ParticleStateInfo &rhs) const {
    if (this->uniqueId_ < rhs.uniqueId_)
      return true;
    else if (this->uniqueId_ > rhs.uniqueId_)
      return false;
    if (this->pid_information_ < rhs.pid_information_)
      return true;
    else if (this->pid_information_ > rhs.pid_information_)
      return false;
    if (this->coherent < rhs.coherent)
      return true;
    else if (this->coherent > rhs.coherent)
      return false;
    if (this->spin_information_ < rhs.spin_information_)
      return true;

    return false;
  }
  bool operator>(const ParticleStateInfo &rhs) const { return (rhs < *this); }

  friend std::ostream &operator<<(std::ostream &os,
                                  const ParticleStateInfo &rhs) {
    os << "unique id: " << rhs.uniqueId_ << std::endl;
    os << "name: " << rhs.pid_information_.name_ << std::endl;
    os << "pid: " << rhs.pid_information_.particleId_ << std::endl;
    os << "J: " << rhs.spin_information_.GetNumerator() << "/"
       << rhs.spin_information_.GetDenominator() << "("
       << rhs.spin_information_.GetZNumerator() << ")";
    if (rhs.coherent)
      os << " coherent" << std::endl;
    return os;
  }
};

} // ns::ComPWA

#endif
