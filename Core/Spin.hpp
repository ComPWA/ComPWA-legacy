// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_SPIN_HPP_
#define COMPWA_SPIN_HPP_

#include <ostream>

#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"
#include "Core/Utils.hpp"

namespace ComPWA {

///
/// Class that represents Fractions limited to denominators of 1 and 2!
///
class Fraction {
public:
  Fraction(int N, unsigned D) : Numerator(N), Denominator(D) { simplify(); }
  Fraction() : Fraction(0, 1) {}
  Fraction(double x) {
    int Signum(1);
    if (x < 0) {
      Signum = -1;
      x = -x;
    }
    unsigned int Num = (unsigned int)x;
    unsigned int Denom(1);
    double AfterDecimal(x - (double)Num);
    if (std::abs(AfterDecimal - 0.5) < 0.01) {
      // its a half-integral spin
      Num = Num * 2 + 1;
      Denom = 2;
    }
    if (Num == 0)
      Signum = 1;
    Numerator = Signum * Num;
    Denominator = Denom;
    simplify();
  }

  // Not the fastest implementation, but speed is irrelevant here
  Fraction &operator-=(const Fraction &rhs) {
    Numerator = Numerator * rhs.Denominator - rhs.Numerator * Denominator;
    Denominator *= rhs.Denominator;
    simplify();
    return *this;
  }

  friend Fraction operator-(Fraction lhs, const Fraction &rhs) {
    lhs -= rhs;
    return lhs;
  }

  operator double() const { return 1.0 * Numerator / Denominator; }

  int getNumerator() const { return Numerator; }
  unsigned int getDenominator() const { return Denominator; }

private:
  void simplify() {
    if (Denominator == 0)
      throw BadParameter(
          "Creating a Fraction with a denominator of 0 is invalid!");
    if (Numerator == 0) {
      Denominator = 1;
      return;
    }
    if (Denominator > 1) {
      unsigned int gcd =
          ComPWA::Utils::greatestCommonDivisor(Numerator, Denominator);
      Numerator /= gcd;
      Denominator /= gcd;
    }
  }

  int Numerator;
  unsigned int Denominator;
};

///
/// \class Spin
/// Spin class for integer and half-integer spins.
///
class Spin {

public:
  // main constructor (all other constructors should call this one)
  Spin(double Magnitude, double Projection) : Spin(Magnitude) {
    auto ProjFraction = Fraction(Projection);
    if (Denominator != ProjFraction.getDenominator()) {
      std::stringstream ss;
      ss << "Spin(): Spin magnitude (" << Numerator << "/" << Denominator
         << ") and projection (" << ProjFraction.getNumerator() << "/"
         << ProjFraction.getDenominator()
         << ") are incompatible (different denominator)";
      throw BadParameter(ss.str());
    }

    if (Numerator < std::abs(ProjFraction.getNumerator()))
      throw BadParameter(
          "Spin(): Spin projection cannot be larger than the magnitude.");

    IsProjectionSet = true;
    ProjectionNumerator = ProjFraction.getNumerator();
  }

  // main constructor (all other constructors should call this one)
  Spin(double Magnitude) : ProjectionNumerator(0), IsProjectionSet(false) {
    if (Magnitude < 0.0)
      throw BadParameter("Spin(): Magnitude is negative!");
    auto MagFraction = Fraction(Magnitude);

    if (MagFraction.getDenominator() > 2)
      throw BadParameter(
          "Spin(): Denominator of Spin has to be either 1 or 2.");

    Numerator = (unsigned int)MagFraction.getNumerator();
    Denominator = MagFraction.getDenominator();
  }

  // Spin 0 is default
  Spin() : Spin(0.0, 0.0){};

  operator unsigned int() const {
    if (!isIntegralSpin())
      throw std::runtime_error(
          "Spin::operator unsigned int(): Invalid cast of half-integral spin "
          "to integral spin!");
    return Numerator / Denominator;
  }
  operator double() const { return 1.0 * Numerator / Denominator; }

  Fraction getMagnitude() const { return Fraction(Numerator, Denominator); }

  Fraction getProjection() const {
    if (!IsProjectionSet)
      throw std::runtime_error("Spin::getProjection(): Projection is not set!");
    return Fraction(ProjectionNumerator, Denominator);
  }

  bool isIntegralSpin() const { return Denominator == 1; }
  bool isProjectionSet() const { return IsProjectionSet; }

  friend std::ostream &operator<<(std::ostream &os, Spin s) {
    os << "J=" << s.Numerator;
    if (s.Denominator == 2)
      os << "/" << s.Denominator;
    os << " (z=" << s.ProjectionNumerator << ")";
    return os;
  }

private:
  unsigned int Numerator;
  unsigned int Denominator;
  int ProjectionNumerator;
  bool IsProjectionSet;
};

} // namespace ComPWA

#endif
