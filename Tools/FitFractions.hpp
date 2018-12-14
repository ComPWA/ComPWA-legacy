// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <memory>
#include <vector>

#include "Physics/Amplitude.hpp"
#include "Physics/CoherentIntensity.hpp"
#include "Tools/Integration.hpp"

#ifndef COMPWA_TOOLS_FITFRACTIONS_HPP_
#define COMPWA_TOOLS_FITFRACTIONS_HPP_

namespace ComPWA {
namespace Tools {

/// Calculates the fit fractions using the formula:
/// \f[
///  f_i = \frac{|c_i|^2 \int A_i A_i^*}{\int \sum c_l c_m^* A_l A_m^*}
/// \f]
/// The \f$c_i\f$ are the complex coefficient of the amplitudes \f$A_i\f$ and
/// the denominator is the integral over the whole intensity.
/// The integrals are performed via Monte Carlo integration method.
ComPWA::ParameterList calculateFitFractions(
    std::shared_ptr<ComPWA::Physics::CoherentIntensity> intensity,
    std::vector<DataPoint> sample) {
  LOG(INFO) << "calculating fit fractions...";
  ComPWA::ParameterList FitFractionsList;

  // calculate denominator
  double IntegralDenominator = ComPWA::Tools::Integral(intensity, sample);

  /*TODO: ultimately we want to have all combinations of amplitudes
  A_i x A_j*

  -so loop over the amplitudes and in the second loop start with the same
  index
  -in that way we will only calculate one half of the matrix, plus the
  diagonal
  -then we need a complexconjugate amplitudedecorator, which we can
  use here to perform the operation!
   */

  // calculate nominators
  for (auto x : intensity->getAmplitudes()) {
    std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>> vec;
    vec.push_back(x);
    auto ci = std::make_shared<ComPWA::Physics::CoherentIntensity>(
        "TempIntensity", vec);
    double IntegralNumerator = ComPWA::Tools::Integral(ci, sample);

    double ffVal = IntegralNumerator / IntegralDenominator;
    LOG(TRACE) << "calculateFitFractions(): fit fraction for (" << x << ") is "
               << ffVal;

    FitFractionsList.addParameter(std::make_shared<ComPWA::FitParameter>(
        std::dynamic_pointer_cast<ComPWA::Physics::NamedAmplitude>(x)
            ->getName(),
        ffVal, 0.0));
  }
  LOG(INFO) << "finished fit fraction calculation!";
  return FitFractionsList;
}

} // namespace Tools
} // namespace ComPWA

#endif
