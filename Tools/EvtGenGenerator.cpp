// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "EvtGenGenerator.hpp"

#include "Core/Properties.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"

#include "ThirdParty/EvtGen/EvtGenKine.hh"
#include "ThirdParty/EvtGen/EvtRandom.hh"
#include "ThirdParty/EvtGen/EvtVector4R.hh"

namespace ComPWA {
namespace Tools {

EvtGenGenerator::EvtGenGenerator(const ComPWA::FourMomentum &CMSP4_,
                                 const std::vector<double> &FinalStateMasses_,
                                 unsigned int seed)
    : CMSP4(CMSP4_), FinalStateMasses(FinalStateMasses_),
      RandomEngine(new EvtGenStdRandomEngine(seed)) {
  if (FinalStateMasses.size() < 2)
    throw std::runtime_error("EvtGenGenerator::EvtGenGenerator() | at least "
                             "two final state particles are required!");
  EvtRandom::setRandomEngine(RandomEngine);
}

EvtGenGenerator::EvtGenGenerator(
    const Physics::ParticleStateTransitionKinematicsInfo &KinematicsInfo,
    unsigned int seed)
    : EvtGenGenerator(KinematicsInfo.getInitialStateFourMomentum(),
                      KinematicsInfo.getFinalStateMasses(), seed) {}

EvtGenGenerator::~EvtGenGenerator() { delete RandomEngine; }

ComPWA::Event EvtGenGenerator::generate() {
  ComPWA::Event evt;

  std::vector<EvtVector4R> FourVectors(FinalStateMasses.size());

  double weight =
      EvtGenKine::PhaseSpace(FinalStateMasses.size(), &FinalStateMasses[0],
                             &FourVectors[0], CMSP4.invMass());
  evt.Weight = weight;

  for (auto const &p4 : FourVectors) {
    evt.ParticleList.push_back(
        Particle(p4.get(1), p4.get(2), p4.get(3), p4.get(0)));
  }
  return evt;
}

void EvtGenGenerator::setSeed(unsigned int seed) {
  RandomEngine->setSeed(seed);
}

unsigned int EvtGenGenerator::getSeed() const {
  return RandomEngine->getSeed();
}

double EvtGenGenerator::uniform(double min, double max) {
  return EvtRandom::Flat(min, max);
}

EvtGenStdRandomEngine::EvtGenStdRandomEngine(unsigned int seed_)
    : MersenneTwisterRandomGenerator(seed_), UniformDistribution(0.0, 1.0),
      seed(seed_) {}

void EvtGenStdRandomEngine::setSeed(unsigned int seed_) {
  seed = seed_;
  MersenneTwisterRandomGenerator.seed(seed);
}
unsigned int EvtGenStdRandomEngine::getSeed() const { return seed; }

double EvtGenStdRandomEngine::random() {
  return UniformDistribution(MersenneTwisterRandomGenerator);
}

} // namespace Tools
} // namespace ComPWA
