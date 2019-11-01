// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/PhspVolume.hpp"
#include "Core/Exceptions.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace ComPWA {
namespace Physics {

class IntegrationSample : public std::vector<double> {
public:
  IntegrationSample(double lower, double upper, std::size_t SampleSize = 100);
  const double BinSize;
};

IntegrationSample::IntegrationSample(double lower, double upper,
                                     std::size_t nvalues)
    : BinSize{(upper - lower) / (nvalues - 1)} {
  this->resize(nvalues);
  double bin_size = BinSize;
  double start_value = lower - BinSize;
  std::generate(begin(), end(),
                [&start_value, &bin_size] { return start_value += bin_size; });
}

std::pair<double, double> SRange(double s, std::vector<double> &masses) {
  std::pair<double, double> s_range = {0., std::sqrt(s)};
  // * Compute s lower
  for (std::size_t i = 0; i < masses.size() - 1; ++i)
    s_range.first += masses[i];
  s_range.first *= s_range.first; // squared
  // * Compute s upper
  s_range.second -= masses.back();
  s_range.second *= s_range.second; // squared
  return s_range;
}

std::pair<double, double> PhspVolume(double s, std::vector<double> &masses,
                                     std::size_t SampleSize) {
  if (masses.size() == 2)
    return std::make_pair(PhspVolumeTwoParticles(s, masses[0], masses[1]), 0.);
  else {
    if (masses.size() > 2) {
      // * Create integration sample
      auto s_range = SRange(s, masses);
      auto sample =
          IntegrationSample(s_range.first, s_range.second, SampleSize + 1);
      // * Create profile vector of phasespace for N - 1 particles
      /// @todo Better to reuse this sample per recursion. Now it is
      /// regenerated for each iteration.
      std::vector<double> previousPhsp(sample.size() - 1);
      auto masses_new = masses;
      masses_new.pop_back();
      double mNsq = masses.back() * masses.back();
      for (std::size_t i = 1; i < sample.size(); ++i) { // don't include s' = 0
        previousPhsp[i - 1] = std::sqrt(KallenFunction(s, sample[i], mNsq));
        previousPhsp[i - 1] *=
            PhspVolume(sample[i], masses_new, SampleSize).first;
        previousPhsp[i - 1] *= sample.BinSize;
      }
      // * Integrate that profile vector
      double volume =
          std::accumulate(previousPhsp.begin(), previousPhsp.end(), 0.) * M_PI /
          s;
      return std::make_pair(volume, 0.);
    } else {
      throw ComPWA::ParameterOutOfBound(
          "Cannot compute a phasespace for only one mass in final state");
    }
  }
}

double PhspVolumeTwoParticles(double s, double m1, double m2) {
  return (2 * M_PI) * std::sqrt(KallenFunction(s, m1 * m1, m2 * m2)) / s;
}

double KallenFunction(double x, double y, double z) {
  double result = x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * z * x;
  if (result < 0.)
    return 0.;
  else
    return result;
}

} // namespace Physics
} // namespace ComPWA
