// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_GENERATOR_HPP_
#define COMPWA_GENERATOR_HPP_

#include "Core/Event.hpp"
#include "Core/Random.hpp"

namespace ComPWA {

///
/// \class PhaseSpaceEventGenerator
/// \brief Interface class for PHSP event generators
///
class PhaseSpaceEventGenerator {
public:
  virtual ~PhaseSpaceEventGenerator() = default;
  virtual EventCollection
  generate(unsigned int NumberOfEvents,
           UniformRealNumberGenerator &RandomGenerator) const = 0;
};

} // namespace ComPWA

#endif
