// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_ESTIMATOR_MINLOGLH_SUMMINLOGLH_HPP_
#define COMPWA_ESTIMATOR_MINLOGLH_SUMMINLOGLH_HPP_

#include <memory>
#include <vector>

#include "Core/FitParameter.hpp"
#include "Core/FunctionTree/FunctionTreeEstimatorWrapper.hpp"
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
  SumMinLogLH(std::vector<std::shared_ptr<MinLogLH>> LogLikelihoods_);

  /// Value of minimum log likelhood function.
  double evaluate();

private:
  std::vector<std::shared_ptr<MinLogLH>> LogLikelihoods;
};

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
createSumMinLogLHEstimatorFunctionTree(
    std::vector<std::shared_ptr<ComPWA::FunctionTree::FunctionTree>>
        LogLikelihoods);

std::tuple<ComPWA::FunctionTree::FunctionTreeEstimatorWrapper, FitParameterList>
createSumMinLogLHFunctionTreeEstimator(
    std::vector<ComPWA::FunctionTree::FunctionTreeEstimatorWrapper> Estimators);

} // namespace Estimator
} // namespace ComPWA
#endif
