// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_ESTIMATOR_MINLOGLH_SUMMINLOGLH_HPP_
#define COMPWA_ESTIMATOR_MINLOGLH_SUMMINLOGLH_HPP_

#include <memory>
#include <vector>

#include "Estimator/Estimator.hpp"

namespace ComPWA {
class FunctionTree;

namespace Estimator {
class MinLogLH;

///
/// \class SumMinLogLH
/// Calculates the combined likelihood of multiple MinLogLH.
///
class SumMinLogLH : public Estimator {
public:
  SumMinLogLH(std::vector<std::shared_ptr<MinLogLH>> LogLikelihoods_);

  /// Value of minimum log likelhood function.
  double evaluate() const;

private:
  std::vector<std::shared_ptr<MinLogLH>> LogLikelihoods;
};

std::shared_ptr<FunctionTree> createSumMinLogLHEstimatorFunctionTree(
    std::vector<std::shared_ptr<FunctionTree>> LogLikelihoods);

} // namespace Estimator
} // namespace ComPWA
#endif
