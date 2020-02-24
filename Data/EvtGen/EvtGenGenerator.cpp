// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "EvtGenGenerator.hpp"

#include "Core/Properties.hpp"
#include "Core/Random.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"

#include "ThirdParty/EvtGen/EvtGenKine.hh"
#include "ThirdParty/EvtGen/EvtRandom.hh"
#include "ThirdParty/EvtGen/EvtVector4R.hh"

namespace ComPWA {
namespace Data {
namespace EvtGen {

EvtGenGenerator::EvtGenGenerator(const ComPWA::FourMomentum &CMSP4_,
                                 const std::vector<double> &FinalStateMasses_)
    : CMSP4(CMSP4_), FinalStateMasses(FinalStateMasses_),
      RandomEngine(new EvtGenStdRandomEngine()) {
  if (FinalStateMasses.size() < 2)
    throw std::runtime_error("EvtGenGenerator::EvtGenGenerator() | at least "
                             "two final state particles are required!");
  EvtRandom::setRandomEngine(RandomEngine.get());
}

EvtGenGenerator::EvtGenGenerator(
    const Physics::ParticleStateTransitionKinematicsInfo &KinematicsInfo)
    : EvtGenGenerator(KinematicsInfo.getInitialStateFourMomentum(),
                      KinematicsInfo.getFinalStateMasses()) {}

ComPWA::Event EvtGenGenerator::generate(UniformRealNumberGenerator &gen) const {
  RandomEngine->setRandomNumberGenerator(gen);
  ComPWA::Event evt;

  std::vector<EvtVector4R> FourVectors(FinalStateMasses.size());

  double weight = EvtGenKine::PhaseSpace(
      FinalStateMasses.size(), (double *)(&FinalStateMasses[0]), // const cast
      &FourVectors[0], CMSP4.invariantMass());
  evt.Weight = weight;

  for (auto const &p4 : FourVectors) {
    evt.FourMomenta.push_back(
        FourMomentum(p4.get(1), p4.get(2), p4.get(3), p4.get(0)));
  }
  return evt;
}

EvtGenStdRandomEngine::EvtGenStdRandomEngine() : NumberGenerator(nullptr) {}

void EvtGenStdRandomEngine::setRandomNumberGenerator(
    UniformRealNumberGenerator &NumberGenerator_) {
  NumberGenerator = &NumberGenerator_;
}
double EvtGenStdRandomEngine::random() { return NumberGenerator->operator()(); }

} // namespace EvtGen
} // namespace Data
} // namespace ComPWA
