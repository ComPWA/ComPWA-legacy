// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Some useful functions for Monte-Carlo event generation.
///

#ifndef COMPWA_DATA_GENERATE_HPP_
#define COMPWA_DATA_GENERATE_HPP_

#include <memory>

#include "Core/Event.hpp"
#include "Core/Function.hpp"

namespace ComPWA {
class UniformRealNumberGenerator;
class Kinematics;
class PhaseSpaceEventGenerator;

namespace Data {

ComPWA::EventList generate(unsigned int NumberOfEvents,
                           const ComPWA::Kinematics &Kinematics,
                           ComPWA::Intensity &Intensity,
                           ComPWA::UniformRealNumberGenerator &RandomGenerator);

ComPWA::EventList generate(unsigned int NumberOfEvents,
                           const ComPWA::Kinematics &Kinematics,
                           ComPWA::UniformRealNumberGenerator &RandomGenerator,
                           ComPWA::Intensity &Intensity,
                           const std::vector<ComPWA::Event> &phsp,
                           const std::vector<ComPWA::Event> &phspTrue);

inline ComPWA::EventList
generate(unsigned int NumberOfEvents, const ComPWA::Kinematics &Kinematics,
         ComPWA::UniformRealNumberGenerator &RandomGenerator,
         ComPWA::Intensity &Intensity, const std::vector<ComPWA::Event> &phsp) {
  return generate(NumberOfEvents, Kinematics, RandomGenerator, Intensity, phsp,
                  phsp);
}

std::vector<Event>
generatePhsp(unsigned int nEvents,
             const ComPWA::PhaseSpaceEventGenerator &Generator,
             ComPWA::UniformRealNumberGenerator &RandomGenerator);

ComPWA::EventList generateImportanceSampledPhsp(
    unsigned int NumberOfEvents, const ComPWA::Kinematics &Kinematics,
    const ComPWA::PhaseSpaceEventGenerator &Generator,
    ComPWA::Intensity &Intensity,
    ComPWA::UniformRealNumberGenerator &RandomGenerator);

} // namespace Data
} // namespace ComPWA
#endif
