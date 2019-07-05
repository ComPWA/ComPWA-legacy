// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Some useful functions for Monte-Carlo event generation.
///

#ifndef COMPWA_TOOLS_GENERATE_HPP_
#define COMPWA_TOOLS_GENERATE_HPP_

#include <memory>

#include "Core/Event.hpp"
#include "Core/Function.hpp"

namespace ComPWA {

class Kinematics;
class Generator;

namespace Tools {

std::vector<ComPWA::Event>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::Intensity> Intensity);

std::vector<ComPWA::Event>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::Intensity> Intensity,
         const std::vector<ComPWA::Event> &phsp,
         const std::vector<ComPWA::Event> &phspTrue);

std::vector<ComPWA::Event>
generate(unsigned int NumberOfEvents,
         std::shared_ptr<ComPWA::Kinematics> Kinematics,
         std::shared_ptr<ComPWA::Generator> Generator,
         std::shared_ptr<ComPWA::Intensity> Intensity,
         const std::vector<ComPWA::Event> &phsp) {
  return generate(NumberOfEvents, Kinematics, Generator, Intensity, phsp, phsp);
}

std::vector<ComPWA::Event> generatePhsp(unsigned int nEvents,
                                        std::shared_ptr<ComPWA::Generator> gen);

std::vector<ComPWA::Event>
generateImportanceSampledPhsp(unsigned int NumberOfEvents,
                              std::shared_ptr<ComPWA::Kinematics> Kinematics,
                              std::shared_ptr<ComPWA::Generator> Generator,
                              std::shared_ptr<ComPWA::Intensity> Intensity);

} // namespace Tools
} // namespace ComPWA
#endif
