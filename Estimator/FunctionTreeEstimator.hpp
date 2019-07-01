// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_ESTIMATOR_FUNCTIONTREEESTIMATOR_HPP_
#define COMPWA_ESTIMATOR_FUNCTIONTREEESTIMATOR_HPP_

#include <memory>

#include "Core/FunctionTree.hpp"
#include "Estimator/Estimator.hpp"

namespace ComPWA {
namespace Estimator {

class FunctionTreeEstimator : public Estimator<double> {
public:
  FunctionTreeEstimator(std::shared_ptr<FunctionTree> functiontree)
      : EvaluationTree(functiontree) {
    if (!EvaluationTree) {
      throw std::runtime_error(
          "FunctionTreeEstimator::FunctionTreeEstimator(): "
          "FunctionTree is empty!");
    }
    EvaluationTree->parameter();
    if (!EvaluationTree->sanityCheck()) {
      throw std::runtime_error(
          "FunctionTreeEstimator::FunctionTreeEstimator(): Tree has structural "
          "problems. Sanity check not passed!");
    }
  }

  double evaluate() final {
    auto EstimatedValue =
        std::dynamic_pointer_cast<Value<double>>(EvaluationTree->parameter());
    return EstimatedValue->value();
  }

  void updateParametersFrom(const std::vector<double> &pars){};

  std::vector<double> getParameters() const;

  std::string print(int level) { return EvaluationTree->head()->print(level); }

private:
  std::shared_ptr<FunctionTree> EvaluationTree;
};

} // namespace Estimator
} // namespace ComPWA

#endif
