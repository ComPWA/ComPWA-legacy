// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_ESTIMATOR_MINLOGLH_SUMMINLOGLH_HPP_
#define COMPWA_ESTIMATOR_MINLOGLH_SUMMINLOGLH_HPP_

#include "Core/FitParameter.hpp"
#include "Core/FunctionTree/FunctionTreeEstimator.hpp"
#include "Estimator/Estimator.hpp"

namespace ComPWA {
namespace FunctionTree {
class FunctionTree;
}

namespace Estimator {
class MinLogLH;

///
/// \class SumMinLogLH
/// Calculates the combined likelihood of multiple MinLogLH.
///
class SumMinLogLH : public Estimator<double> {
public:
  SumMinLogLH(std::vector<std::shared_ptr<Estimator>> Estimators);

  /// Value of minimum log likelihood function.
  double evaluate() noexcept final;

  void updateParametersFrom(const std::vector<double> &params) final;

  std::vector<ComPWA::Parameter> getParameters() const final;

private:
  std::vector<std::shared_ptr<Estimator>> LogLikelihoods;
};

std::tuple<ComPWA::FunctionTree::FunctionTreeEstimator, FitParameterList>
createSumMinLogLHFunctionTreeEstimator(
    std::vector<std::pair<ComPWA::FunctionTree::FunctionTreeEstimator,
                          FitParameterList>>
        Estimators);

} // namespace Estimator
} // namespace ComPWA
#endif
