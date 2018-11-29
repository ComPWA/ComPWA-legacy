/*
 * EvtGenGenerator.cpp
 *
 *  Created on: Nov 21, 2018
 *      Author: steve
 */

#include "EvtGenGenerator.hpp"

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
  EvtRandom::setRandomEngine(RandomEngine);
}

EvtGenGenerator::EvtGenGenerator(std::shared_ptr<PartList> partL,
                                 std::shared_ptr<Kinematics> kin,
                                 unsigned int seed)
    : RandomEngine(new EvtGenStdRandomEngine(seed)) {
  EvtRandom::setRandomEngine(RandomEngine);
  auto const &KinProps(kin->getKinematicsProperties());
  auto finalS = KinProps.FinalState;
  auto initialS = KinProps.InitialState;
  unsigned int nPart = finalS.size();
  if (nPart < 2)
    throw std::runtime_error(
        "EvtGenGenerator::RootGenerator() | one particle is not enough!");
  CMSP4 = KinProps.InitialStateP4;
  for (auto ParticlePid : finalS) { // particle 0 is mother particle
    FinalStateMasses.push_back(FindParticle(partL, ParticlePid).GetMass());
  }
}

EvtGenGenerator::~EvtGenGenerator() { delete RandomEngine; }

ComPWA::Event EvtGenGenerator::generate() {
  ComPWA::Event evt;

  std::vector<EvtVector4R> FourVectors(FinalStateMasses.size());

  double weight =
      EvtGenKine::PhaseSpace(FinalStateMasses.size(), &FinalStateMasses[0],
                             &FourVectors[0], CMSP4.invMass());
  evt.setWeight(weight);

  // final boost of all particles
  for (auto const &p4 : FourVectors) {
    evt.addParticle(Particle(p4.get(1), p4.get(2), p4.get(3), p4.get(0)));
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
