// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Estimator/MinLogLH/SumMinLogLH.hpp"
#include "Core/Event.hpp"
#include "Core/FitResult.hpp"
#include "Core/FunctionTree/FunctionTree.hpp"
#include "Core/FunctionTree/ParameterList.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Particle.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"

namespace ComPWA {
namespace Estimator {

using namespace ComPWA::FunctionTree;

SumMinLogLH::SumMinLogLH(std::vector<std::shared_ptr<Estimator>> Estimators)
    : LogLikelihoods(Estimators) {}

double SumMinLogLH::evaluate() noexcept {
  double lh(0.0);
  for (auto x : LogLikelihoods)
    lh += x->evaluate();
  return lh;
}

void SumMinLogLH::updateParametersFrom(const std::vector<double> &params) {
  auto CurrentSubRangeBegin = params.begin();
  for (auto LLH : LogLikelihoods) {
    size_t SubRangeSize = LLH->getParameters().size();
    std::vector<double> SubParameterRange;
    SubParameterRange.reserve(SubRangeSize);
    std::copy(CurrentSubRangeBegin, CurrentSubRangeBegin + SubRangeSize,
              std::back_inserter(SubParameterRange));
    CurrentSubRangeBegin += SubRangeSize;
    LLH->updateParametersFrom(SubParameterRange);
  }
}

std::vector<ComPWA::Parameter> SumMinLogLH::getParameters() const {
  std::vector<ComPWA::Parameter> Parameters;
  for (auto x : LogLikelihoods) {
    auto pars = x->getParameters();
    Parameters.insert(Parameters.end(), pars.begin(), pars.end());
  }
  return Parameters;
}

// This function is one part, which suggest a bad design
std::tuple<FunctionTreeEstimator, FitParameterList>
createSumMinLogLHFunctionTreeEstimator(
    std::vector<std::pair<ComPWA::FunctionTree::FunctionTreeEstimator,
                          FitParameterList>>
        Estimators) {

  unsigned int counter(1);
  FitParameterList TempParameters;
  FitParameterList Pars;
  ParameterList ParList;

  auto EvaluationTree = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      "SumLogLh", std::make_shared<Value<double>>(),
      std::make_shared<AddAll>(ParType::DOUBLE));

  for (auto &x : Estimators) {
    try {
      // we need to change the names of the log likelihoods so that the
      // function tree will be constructed correctly
      x.first.getFunctionTree()->Head->setName("LH_" + std::to_string(counter));
      EvaluationTree->insertTree(x.first.getFunctionTree(), "SumLogLh");
    } catch (std::exception &ex) {
      LOG(ERROR) << "createSumMinLogLHEstimatorFunctionTree(): Construction of "
                    "one or more sub trees has failed! Error: "
                 << ex.what();
    }
    ++counter;
    TempParameters.insert(TempParameters.begin(), x.second.begin(),
                          x.second.end());
  }
  EvaluationTree->Head->fillParameters(ParList);

  for (auto x : ParList.doubleParameters()) {
    auto result = std::find_if(TempParameters.begin(), TempParameters.end(),
                               [x](const ComPWA::FitParameter<double> &fp) {
                                 return fp.Name == x->name();
                               });
    if (TempParameters.end() != result) {
      Pars.push_back(*result);
    }
  }

  EvaluationTree->parameter();

  return std::make_tuple(FunctionTreeEstimator(EvaluationTree, ParList), Pars);
} // namespace Estimator

} // namespace Estimator
} // namespace ComPWA
