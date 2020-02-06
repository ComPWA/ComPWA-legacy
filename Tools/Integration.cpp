// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Integration.hpp"

#include "Core/Logging.hpp"
#include "Data/DataSet.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "ThirdParty/parallelstl/include/pstl/algorithm"
#include "ThirdParty/parallelstl/include/pstl/execution"

namespace ComPWA {
namespace Tools {

struct KahanSummation {
  double sum;
  double correction;
  operator double() const { return sum; }
};

/// KahanSummation keeps track of lost bits and reduced the uncertainty in the
/// summation of many large/small numbers.
/// See https://en.wikipedia.org/wiki/Kahan_summation_algorithm
KahanSummation KahanSum(KahanSummation accumulation, double value) {
  KahanSummation result;
  double y = value - accumulation.correction;
  double t = accumulation.sum + y;
  result.correction = (t - accumulation.sum) - y;
  result.sum = t;
  return result;
}

std::pair<double, double>
integrateWithError(ComPWA::Intensity &intensity,
                   const ComPWA::Data::DataSet &phspsample, double phspVolume) {
  auto WeightSum =
      std::accumulate(phspsample.Weights.begin(), phspsample.Weights.end(),
                      KahanSummation{0., 0.}, KahanSum);

  std::vector<double> Intensities = intensity.evaluate(phspsample.Data);
  std::transform(
      pstl::execution::par_unseq, Intensities.begin(), Intensities.end(),
      phspsample.Weights.begin(), Intensities.begin(),
      [](double intensity, double weight) { return intensity * weight; });

  auto IntensitySum = std::accumulate(Intensities.begin(), Intensities.end(),
                                      KahanSummation{0., 0.}, KahanSum);
  double AvgInt = IntensitySum / WeightSum;
  double Integral = AvgInt * phspVolume;

  // We reuse the intensities vector to store residuals of intensities
  std::transform(pstl::execution::par_unseq, Intensities.begin(),
                 Intensities.end(), Intensities.begin(),
                 [&AvgInt](double intensity) {
                   return (intensity - AvgInt) * (intensity - AvgInt);
                 });
  auto IntensityResidualsSum = std::accumulate(
      Intensities.begin(), Intensities.end(), KahanSummation{0., 0.}, KahanSum);
  double AvgIntResSq = IntensityResidualsSum / (WeightSum - 1);
  double IntegralErrorSq = AvgIntResSq * phspVolume * phspVolume / WeightSum;

  return std::make_pair(Integral, std::sqrt(IntegralErrorSq));
}

double integrate(ComPWA::Intensity &intensity,
                 const ComPWA::Data::DataSet &phspsample, double phspVolume) {
  return integrateWithError(intensity, phspsample, phspVolume).first;
}

double maximum(ComPWA::Intensity &intensity,
               const ComPWA::Data::DataSet &sample) {
  if (!sample.Weights.size()) {
    LOG(DEBUG) << "Tools::Maximum(): Maximum can not be determined since "
                  "sample is empty.";
    return 1.0;
  }

  std::vector<double> Intensities = intensity.evaluate(sample.Data);

  std::transform(Intensities.begin(), Intensities.end(), sample.Weights.begin(),
                 Intensities.begin(), [](double Intensity, double Weight) {
                   return Intensity * Weight;
                 });
  // determine maximum
  double max(*std::max_element(Intensities.begin(), Intensities.end()));
  LOG(INFO) << "Tools::Maximum(): found maximum value of " << max;
  return max;
}

} // namespace Tools
} // namespace ComPWA
