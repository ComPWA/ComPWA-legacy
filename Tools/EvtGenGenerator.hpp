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

namespace Physics {
class ParticleStateTransitionKinematicsInfo;
}

namespace Tools {

class EvtGenStdRandomEngine : public EvtRandomEngine {
  std::mt19937 MersenneTwisterRandomGenerator;
  std::uniform_real_distribution<double> UniformDistribution;
  unsigned int seed;

public:
  EvtGenStdRandomEngine(unsigned int seed_);

  void setSeed(unsigned int seed_);
  unsigned int getSeed() const;

  double random();
};

class EvtGenGenerator : public Generator {
  ComPWA::FourMomentum CMSP4;
  std::vector<double> FinalStateMasses;
  ComPWA::Tools::EvtGenStdRandomEngine *RandomEngine;

public:
  EvtGenGenerator(const ComPWA::FourMomentum &CMSP4_,
                  const std::vector<double> &FinalStateMasses_,
                  unsigned int seed);

  /// Constructor: Information on the decay is obtained from Kinematics
  EvtGenGenerator(
      const Physics::ParticleStateTransitionKinematicsInfo &KinematicsInfo,
      unsigned int seed = -1);

  virtual ~EvtGenGenerator();

  ComPWA::Event generate();

  void setSeed(unsigned int seed);

  unsigned int getSeed() const;

  double uniform(double min, double max);

  double gauss(double mu, double sigma) const { return 0; }
};

} // namespace Tools
} // namespace ComPWA

#endif
