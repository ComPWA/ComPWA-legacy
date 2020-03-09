// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef TOOLS_ROOTGENERATOR_HPP_
#define TOOLS_ROOTGENERATOR_HPP_

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "TRandom3.h"

#include "Core/Generator.hpp"
#include "Core/Properties.hpp"

namespace ComPWA {
class UniformRealNumberGenerator;
namespace Physics {
class ParticleStateTransitionKinematicsInfo;
}

namespace Data {
namespace Root {

class RootUniformRealGenerator : public UniformRealNumberGenerator {
  TRandom3 RandomGenerator;
  int Seed;

public:
  RootUniformRealGenerator(int seed = 123456);

  double operator()() final;
  int getSeed() const final;
  void setSeed(int seed) final;
};

class RootGenerator : public PhaseSpaceEventGenerator {
public:
  /// Constructor for a three particle decay with given masses
  RootGenerator(const ComPWA::FourMomentum &CMSP4_,
                const std::vector<double> &FinalStateMasses_,
                const std::vector<ComPWA::pid> &FinalStatePIDs_);

  /// Constructor: Information on the decay is obtained from Kinematics
  RootGenerator(
      const Physics::ParticleStateTransitionKinematicsInfo &KinematicsInfo);

  /// Constructor: Information on the decay is provides via lists of initial and
  /// final states
  RootGenerator(const ComPWA::ParticleList &PartL, std::vector<pid> FinalS,
                std::vector<pid> InitialS);

  virtual ~RootGenerator(){};

  ComPWA::EventCollection
  generate(unsigned int NumberOfEvents,
           UniformRealNumberGenerator &RandomGenerator) const final;

private:
  void init();
  /// These functions are copied from ROOT
  double PDK(double a, double b, double c) const;
  void BoostAlongY(TLorentzVector &vec, double beta_squared) const;

  ComPWA::FourMomentum CMSP4;
  std::vector<double> FinalStateMasses;
  std::vector<ComPWA::pid> FinalStatePIDs;
  double MaximumWeight;
  TVector3 CMSBoostVector;
  // total energy in C.M. minus the sum of the masses
  double CMSEnergyMinusMasses;
};

} // namespace Root
} // namespace Data
} // namespace ComPWA

#endif
