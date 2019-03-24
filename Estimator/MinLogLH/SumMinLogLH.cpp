// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Estimator/MinLogLH/SumMinLogLH.hpp"
#include "Core/Event.hpp"
#include "Core/FitResult.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Kinematics.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Particle.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"

namespace ComPWA {
namespace Estimator {

SumMinLogLH::SumMinLogLH(std::vector<std::shared_ptr<MinLogLH>> LogLikelihoods_)
    : LogLikelihoods(LogLikelihoods_) {}

double SumMinLogLH::evaluate() const {
  double lh(0.0);
  for (auto const x : LogLikelihoods)
    lh += x->evaluate();
  return lh;
}

std::shared_ptr<FunctionTree> createSumMinLogLHEstimatorFunctionTree(
    std::vector<std::shared_ptr<FunctionTree>> LogLikelihoods) {
  auto EvaluationTree = std::make_shared<FunctionTree>(
      "SumLogLh", std::make_shared<Value<double>>(),
      std::make_shared<AddAll>(ParType::DOUBLE));
  unsigned int counter(1);
  for (auto ll : LogLikelihoods) {
    try {
      // we need to change the names of the log likelihoods so that the
      // function tree will be constructed correctly
      ll->head()->setName("LH_" + std::to_string(counter));
      EvaluationTree->insertTree(ll, "SumLogLh");
    } catch (std::exception &ex) {
      LOG(ERROR) << "createSumMinLogLHEstimatorFunctionTree(): Construction of "
                    "one or more sub trees has failed! Error: "
                 << ex.what();
    }
    ++counter;
  }

  EvaluationTree->parameter();
  if (!EvaluationTree->sanityCheck()) {
    throw std::runtime_error(
        "createSumMinLogLHEstimatorFunctionTree(): tree has structural "
        "problems. Sanity check not passed!");
  }
  return EvaluationTree;
}

} // namespace Estimator
} // namespace ComPWA
