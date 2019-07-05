// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "NonResonant.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

std::complex<double> NonResonant::evaluate(const ComPWA::DataPoint &point,
                                           unsigned int pos) const {
  return std::complex<double>(1.0, 0.0);
}

void NonResonant::updateParametersFrom(
    const ComPWA::FunctionTree::ParameterList &list) {}
void NonResonant::addUniqueParametersTo(
    ComPWA::FunctionTree::ParameterList &list) {}
void NonResonant::addFitParametersTo(std::vector<double> &FitParameters) {}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
NonResonant::createFunctionTree(
    const ComPWA::FunctionTree::ParameterList &DataSample, unsigned int pos,
    const std::string &suffix) const {

  auto unitVec = std::make_shared<FunctionTree::Value<std::complex<double>>>(
      std::complex<double>(1, 0));

  return std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      "NonResonant" + suffix, unitVec);
}

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA
