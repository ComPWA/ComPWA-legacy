// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef TOOLS_EVTGENGENERATOR_HPP_
#define TOOLS_EVTGENGENERATOR_HPP_

#include <random>

#include "Core/Event.hpp"
#include "Core/Generator.hpp"
#include "ThirdParty/EvtGen/EvtRandomEngine.hh"

namespace ComPWA {
class UniformRealNumberGenerator;
namespace Physics {
class ParticleStateTransitionKinematicsInfo;
}

namespace Data {
namespace EvtGen {

class EvtGenStdRandomEngine : public EvtRandomEngine {
  // ownership is not taken of this random number generator
  UniformRealNumberGenerator *NumberGenerator;

public:
  EvtGenStdRandomEngine();
  void setRandomNumberGenerator(UniformRealNumberGenerator &NumberGenerator_);

  double random();
};

class EvtGenGenerator : public PhaseSpaceEventGenerator {
  ComPWA::FourMomentum CMSP4;
  std::vector<double> FinalStateMasses;
  std::vector<ComPWA::pid> FinalStatePIDs;
  std::unique_ptr<ComPWA::Data::EvtGen::EvtGenStdRandomEngine> RandomEngine;

public:
  EvtGenGenerator(const ComPWA::FourMomentum &CMSP4_,
                  const std::vector<double> &FinalStateMasses_,
                  const std::vector<ComPWA::pid> &FinalStatePIDs_);

  /// Constructor: Information on the decay is obtained from Kinematics
  EvtGenGenerator(
      const Physics::ParticleStateTransitionKinematicsInfo &KinematicsInfo);

  ComPWA::EventCollection
  generate(unsigned int NumberOfEvents,
           UniformRealNumberGenerator &RandomGenerator) const final;
};

} // namespace EvtGen
} // namespace Data
} // namespace ComPWA

#endif
