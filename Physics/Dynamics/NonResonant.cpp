// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "NonResonant.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
NonResonant::createFunctionTree(
    const ComPWA::FunctionTree::ParameterList &DataSample, unsigned int pos,
    std::string suffix) {

  auto unitVec = std::make_shared<FunctionTree::Value<std::complex<double>>>(
      std::complex<double>(1, 0));

  return std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      "NonResonant" + suffix, unitVec);
}

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA
