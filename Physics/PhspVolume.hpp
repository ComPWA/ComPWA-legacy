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

#include "Core/Random.hpp"
#include "Data/DataSet.hpp"
#include "Tools/Integration.hpp"
#include <vector>

namespace ComPWA {
namespace Physics {

/// Compute phasespace volume of momentum space for an arbitrary number of
/// particles in the final state *using Monte Carlo integration*.
std::pair<double, double> phspVolumeMC(double s, std::vector<double> FSMasses,
                                       size_t NumEvaluations = 100000);

/// Compute phasespace volume of momentum space for an arbitrary number of
/// particles in the final state *using Riemann integration*.
/// @TODO Algorithm might be improved with [Simpson's
/// rule](https://en.wikipedia.org/wiki/Simpson%27s_rule), because we integrate
/// over a function that is polynomial in the limit \f$m_i\rightarrow 0\f$
std::pair<double, double> phspVolumeRiemann(double s,
                                            std::vector<double> &FSMasses,
                                            size_t NumEvaluations = 1000);

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

class PhaseSpaceVolume : public Intensity {
public:
  PhaseSpaceVolume(double ISMass_, std::vector<double> FSMasses_)
      : ISMass(ISMass_), FSMassesSquared(FSMasses_) {
    for (auto &x : FSMassesSquared)
      x = x * x;
  }

  std::vector<double>
  evaluate(const std::vector<std::vector<double>> &points) noexcept;

  void updateParametersFrom(const std::vector<double> &){};
  std::vector<Parameter> getParameters() const { return {}; };

private:
  double ISMass;
  std::vector<double> FSMassesSquared;
};

} // namespace Physics
} // namespace ComPWA

#endif
