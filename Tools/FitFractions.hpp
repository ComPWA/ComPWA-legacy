// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_TOOLS_FITFRACTIONS_HPP_
#define COMPWA_TOOLS_FITFRACTIONS_HPP_

#include "Core/Function.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ComPWA {
struct FitResult;
class Kinematics;
namespace Data {
struct DataSet;
}
namespace Tools {

struct FitFraction {
  std::string Name;
  double Value;
  double Error;
};

using FitFractionList = std::vector<FitFraction>;

using IntensityComponent = std::pair<std::string, std::shared_ptr<Intensity>>;

std::ostream &operator<<(std::ostream &os, const FitFractionList &FFList);

class FitFractions {
public:
  /// Calculates the fit fractions with errors via error propagation from the
  /// covariance matrix. The gradients are calculated via numerical
  /// differentiation:
  /// \f[
  /// f´(x) = \frac{f(x+h) - f(x-h)}{2h} + O(h^2)
  /// \f]
  FitFractionList calculateFitFractionsWithCovarianceErrorPropagation(
      const std::vector<std::pair<IntensityComponent, IntensityComponent>>
          &Components,
      const ComPWA::Data::DataSet &PhspSample, const ComPWA::FitResult &Result);

private:
  // internal helper structure, hidden from user
  struct DerivativeData {
    std::string ParameterName;
    double ValueAtParameterPlusEpsilon;
    double ValueAtParameterMinusEpsilon;
    double StepSize;
  };

  std::map<std::string, std::tuple<double, std::vector<DerivativeData>>>
      IntensityGradientDataMapping;

  std::tuple<std::vector<double>, std::vector<std::vector<double>>>
  buildJacobiAndCovariance(const std::tuple<double, std::vector<DerivativeData>>
                               &NominatorDerivatives,
                           const std::tuple<double, std::vector<DerivativeData>>
                               &DenominatorDerivatives,
                           const ComPWA::FitResult &Result);

  std::tuple<std::string, double, std::vector<FitFractions::DerivativeData>>
  getIntegralData(IntensityComponent IntensComponent,
                  const ComPWA::Data::DataSet &PhspSample,
                  const ComPWA::FitResult &Result);

  std::tuple<double, std::vector<FitFractions::DerivativeData>>
  calculateIntensityIntegralData(ComPWA::Intensity &Intens,
                                 const ComPWA::Data::DataSet &PhspSample,
                                 const ComPWA::FitResult &Result);
};

} // namespace Tools
} // namespace ComPWA

#endif
