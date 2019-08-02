// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_DYNAMICS_NONRESONANT_HPP_
#define COMPWA_PHYSICS_DYNAMICS_NONRESONANT_HPP_

#include "Core/FunctionTree/FunctionTree.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

namespace NonResonant {
std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
createFunctionTree(const ComPWA::FunctionTree::ParameterList &DataSample,
                   unsigned int pos, std::string suffix);
}

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
