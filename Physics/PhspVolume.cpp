// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/PhspVolume.hpp"
#include "Core/Exceptions.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

namespace ComPWA {
namespace Physics {

double PhspVolume(double s, std::vector<double> &masses, size_t nsteps) {
  if (masses.size() == 2)
    return PhspTwoParticles(s, masses[0], masses[1]);
  else {
    if (masses.size() > 2) {
      // * Create integration sample
      /// @TODO: Results in many integration samples when computing 5 particles
      auto s_range = SRange(s, masses);
      auto sample =
          IntegrationSample(s_range.first, s_range.second, nsteps + 1);
      // * Create profile vector of phasespace for N - 1 particles
      std::vector<double> previousPhsp(sample.size() - 1);
      auto masses_new = masses;
      masses_new.pop_back();
      double mNsq = masses.back() * masses.back();
      for (size_t i = 1; i < sample.size(); ++i) { // don't include s' = 0
        previousPhsp[i - 1] = SqrtKallenFunction(s, sample[i], mNsq);
        previousPhsp[i - 1] *= PhspVolume(sample[i], masses_new, nsteps);
        previousPhsp[i - 1] *= sample.BinSize;
      }
      // * Integrate that profile vector
      return std::accumulate(previousPhsp.begin(), previousPhsp.end(), 0.) *
             M_PI / s;
    } else {
      throw ComPWA::ParameterOutOfBound(
          "Cannot compute a phasespace for only one mass in final state");
      return 0.;
    }
  }
}

double PhspTwoParticles(double s, double m1, double m2) {
  return (2 * M_PI) * ThresholdFunction(s, m1, m2);
}

double ThresholdFunction(double s, double m1, double m2) {
  return SqrtKallenFunction(s, m1 * m1, m2 * m2) / s;
}

double SqrtKallenFunction(double x, double y, double z) {
  double result = KallenFunction(x, y, z);
  if (result < 0.)
    return 0.;
  else
    return std::sqrt(result);
}

double KallenFunction(double x, double y, double z) {
  return x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * z * x;
}

std::pair<double, double> SRange(double s, std::vector<double> &masses) {
  std::pair<double, double> s_range = {0., std::sqrt(s)};
  // * Compute s lower
  for (size_t i = 0; i < masses.size() - 1; ++i)
    s_range.first += masses[i];
  s_range.first *= s_range.first; // squared
  // * Compute s upper
  s_range.second -= masses.back();
  s_range.second *= s_range.second; // squared
  return s_range;
}

/// Generate a vector of size \p nvalues of equally spaced values lying between
/// and including \p lower and \p upper. The bin size is stored as a data
/// member, because it is used repeatedly when computing a Rieman integral.
/// @TODO: Move to general tools and make into a template. Note that you'll run
/// into trouble at `generate`.
IntegrationSample::IntegrationSample(double lower, double upper, size_t nvalues)
    : BinSize{(upper - lower) / (nvalues - 1)} {
  this->resize(nvalues);
  double bin_size = BinSize;
  double start_value = lower - BinSize;
  std::generate(begin(), end(),
                [&start_value, &bin_size] { return start_value += bin_size; });
}

} // namespace Physics
} // namespace ComPWA
