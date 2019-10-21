// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

/// \file
/// Function to calculate the available phase space volume of arbitrary decays.
/// The algorithm follows the recursive integral computation outlined in [this
/// paper](http://theory.gsi.de/~knoll/Lecture-notes/1-kinematic.pdf), see
/// pp.&nbsp;6–7.

#ifndef PhspVolume_h
#define PhspVolume_h

#include <vector>

namespace ComPWA {
namespace Physics {

class IntegrationSample : public std::vector<double> {
public:
  IntegrationSample(double lower, double upper, size_t nsteps = 100);
  const double BinSize;
};

/// Original [Källén
/// function](https://en.wikipedia.org/wiki/K%C3%A4ll%C3%A9n_function), that is,
/// not having square values in its argument. We use this function instead of
/// the one that can be factorised (see [Heron's
/// formula](https://en.wikipedia.org/wiki/Heron%27s_formula)), because we need
/// to enter \f$s\f$ without taking its square root.
double KallenFunction(double x, double y, double z);
double SqrtKallenFunction(double x, double y, double z);

/// Threshold function. Can be interpreted as the particle velocity at a given
/// center-of-mass energy \p s in case of equal masses (\p m1 = \p m2).
/// Vanishes below threshold.
double ThresholdFunction(double s, double m1, double m2);

/// Phase space element for a two particle decay. An analytic solution exists
/// only for the volume of the phasespace of two-particle decays.
double PhspTwoParticles(double s, double m1, double m2);

std::pair<double, double> SRange(double s, std::vector<double> &masses);

double PhspVolume(double s, std::vector<double> &masses, size_t nsteps = 100);

} // namespace Physics
} // namespace ComPWA

#endif
