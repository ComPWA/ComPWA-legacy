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

double integrate(ComPWA::Intensity &intensity,
                 const ComPWA::Data::DataSet &phspsample, double phspVolume) {

  std::vector<double> Intensities = intensity.evaluate(phspsample.Data);
  std::transform(
      pstl::execution::par_unseq, Intensities.begin(), Intensities.end(),
      phspsample.Weights.begin(), Intensities.begin(),
      [](double intensity, double weight) { return intensity * weight; });
  double IntensitySum(
      std::accumulate(Intensities.begin(), Intensities.end(), 0.0));
  double WeightSum(std::accumulate(phspsample.Weights.begin(),
                                   phspsample.Weights.end(), 0.0));

  return (IntensitySum * phspVolume / WeightSum);
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
