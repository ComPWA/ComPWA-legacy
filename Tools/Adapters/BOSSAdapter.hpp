// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef TOOLS_ADAPTERS_BOSSADAPTER_HPP_
#define TOOLS_ADAPTERS_BOSSADAPTER_HPP_

#include <memory>
#include <utility>
#include <vector>

#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

namespace ComPWA {
namespace Tools {
namespace Adapter {
namespace BOSS {

std::pair<FunctionTree::FunctionTreeIntensity,
          Physics::HelicityFormalism::HelicityKinematics>
createHelicityModel(const char *modelXMLFile, int seed,
                    const std::vector<pid> &initialState,
                    const std::vector<pid> &finalState,
                    const char *particleXMLFile);
}
} // namespace Adapter
} // namespace Tools
} // namespace ComPWA

#endif
