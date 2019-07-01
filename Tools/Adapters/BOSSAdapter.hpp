// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef TOOLS_ADAPTERS_BOSSADAPTER_HPP_
#define TOOLS_ADAPTERS_BOSSADAPTER_HPP_

#include <memory>
#include <utility>
#include <vector>

#include "Core/Function.hpp"

namespace ComPWA {

class Kinematics;

namespace Tools {
namespace Adapter {
namespace BOSS {

std::pair<std::shared_ptr<Intensity>, std::shared_ptr<Kinematics>>
createHelicityModel(const char *modelXMLFile, int seed,
                    const std::vector<int> &initialState,
                    const std::vector<int> &finalState,
                    const char *particleXMLFile);
}
} // namespace Adapter
} // namespace Tools
} // namespace ComPWA

#endif
