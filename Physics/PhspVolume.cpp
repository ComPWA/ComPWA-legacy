// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/PhspVolume.hpp"
#include "Core/Exceptions.hpp"
#include "Tools/Integration.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

namespace ComPWA {
namespace Physics {

std::pair<double, double>
phspVolumeRiemann(double s, std::vector<double> &masses, size_t nsteps) {
  if (masses.size() == 2)
    return std::make_pair(PhspTwoParticles(s, masses[0], masses[1]), 0.);
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
        previousPhsp[i - 1] *=
            phspVolumeRiemann(sample[i], masses_new, nsteps).first;
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

std::pair<double, double> phspVolumeMC(double s, std::vector<double> FSMasses,
                                       size_t NumEvaluations) {
  PhaseSpaceVolume PSVolumeFunction(s, FSMasses);

  // define the intervals (identical to the definition above)
  std::vector<std::pair<double, double>> SIntervals;
  auto FSMassBegin = FSMasses.begin() + 2;
  auto FSMassEnd = FSMasses.end();
  double upper(s);
  double TempMassSum(std::accumulate(FSMasses.begin(), FSMasses.end(), 0.0));
  for (auto x = FSMassBegin; x != FSMassEnd; ++x) {
    upper = std::pow(std::sqrt(upper) - *x, 2);
    TempMassSum -= *x;
    SIntervals.push_back(std::make_pair(std::pow(TempMassSum, 2), upper));
  }
  std::reverse(std::begin(SIntervals), std::end(SIntervals));

  // generate "points" uniform in the s intervals
  ComPWA::Data::DataSet Data;
  StdUniformRealGenerator Random(1234);
  std::vector<double> point(SIntervals.size());
  for (size_t i = 0; i < NumEvaluations; ++i) {
    std::transform(SIntervals.begin(), SIntervals.end(), point.begin(),
                   [&Random](const std::pair<double, double> &limits) {
                     double r = Random();
                     return limits.first + r * (limits.second - limits.first);
                   });
    Data.Data.push_back(point);
    Data.Weights.push_back(1.0);
  }

  double HyperCubeVolume(1.0);
  for (auto const &x : SIntervals) {
    HyperCubeVolume *= (x.second - x.first);
  }

  // then just call the numerical integration
  return ComPWA::Tools::integrateWithError(PSVolumeFunction, Data,
                                           HyperCubeVolume);
}

std::vector<double> PhaseSpaceVolume::evaluate(
    const std::vector<std::vector<double>> &points) noexcept {
  std::vector<double> results;
  for (auto point : points) {
    double result(0.0);
    if (point.size() > 0) {
      result = KallenFunction(ISMass, point.back(), FSMassesSquared.back());
      for (size_t i = 1; i < point.size(); ++i) {
        if (point[i] < point[i - 1]) {
          result = 0.0;
          break;
        }
        result *=
            KallenFunction(point[i], point[i - 1], FSMassesSquared[i + 1]);
      }
      result *=
          KallenFunction(point.front(), FSMassesSquared[0], FSMassesSquared[1]);
      if (result < 0) {
        result = 0.0;
      }
    } else {
      result = KallenFunction(ISMass, FSMassesSquared[0], FSMassesSquared[1]);
    }

    if (result > 0.0) {
      result = 2.0 * M_PI / ISMass * std::sqrt(result);
      for (auto const &x : point) {
        result *= M_PI / x;
      }
    }

    results.push_back(result);
  }
  return results;
}

} // namespace Physics
} // namespace ComPWA
