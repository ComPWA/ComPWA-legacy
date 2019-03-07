// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Some functions to quickly create an intensity object
///

#ifndef COMPWA_TOOLS_CREATEINTENSITY_HPP_
#define COMPWA_TOOLS_CREATEINTENSITY_HPP_

#include <memory>
#include "Core/IIntensity.hpp"
#include "Core/IKinematics.hpp"

namespace ComPWA {
namespace Tools {

double intensityMaximum(std::shared_ptr<ComPWA::IIntensity> amp, int size);

/// Create Intensity and Kinematics class
///
/// This function is used to interface with ComPWA from ThirdParty code.
/// This idea is that is interface in backward compatible to old compilers and
/// old stdlib (prior C++11). Due to an ABI change in stdlib we use char*
/// instead of std::strings.
std::pair<std::shared_ptr<ComPWA::IIntensity>, std::shared_ptr<IKinematics>>
createHelicityModel(const char* modelFile, int seed, int mcPrecision = 1000000,
                    const char* logLv = "FATAL");

}  // namespace Tools
}  // namespace ComPWA
#endif
