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

ComPWA::EventCollection
generate(unsigned int NumberOfEvents, const ComPWA::Kinematics &Kinematics,
         const ComPWA::PhaseSpaceEventGenerator &Generator,
         ComPWA::Intensity &Intensity,
         ComPWA::UniformRealNumberGenerator &RandomGenerator);

ComPWA::EventCollection
generate(unsigned int NumberOfEvents, const ComPWA::Kinematics &Kinematics,
         ComPWA::UniformRealNumberGenerator &RandomGenerator,
         ComPWA::Intensity &Intensity, const EventCollection &PhspSample,
         const EventCollection &PhspSampleTrue);

inline ComPWA::EventCollection
generate(unsigned int NumberOfEvents, const ComPWA::Kinematics &Kinematics,
         ComPWA::UniformRealNumberGenerator &RandomGenerator,
         ComPWA::Intensity &Intensity, const EventCollection &PhspSample) {
  return generate(NumberOfEvents, Kinematics, RandomGenerator, Intensity,
                  PhspSample, PhspSample);
}

ComPWA::EventCollection
generatePhsp(unsigned int NumberOfEvents,
             const ComPWA::PhaseSpaceEventGenerator &Generator,
             ComPWA::UniformRealNumberGenerator &RandomGenerator);

ComPWA::EventCollection generateImportanceSampledPhsp(
    unsigned int NumberOfEvents, const ComPWA::Kinematics &Kinematics,
    const ComPWA::PhaseSpaceEventGenerator &Generator,
    ComPWA::Intensity &Intensity,
    ComPWA::UniformRealNumberGenerator &RandomGenerator);

} // namespace Data
} // namespace ComPWA
#endif
