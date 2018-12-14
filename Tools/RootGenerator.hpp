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

namespace Physics {
class ParticleStateTransitionKinematicsInfo;
}

namespace Tools {

class RootGenerator : public Generator {
  /// These functions are copied from ROOT
  double PDK(double a, double b, double c) const;
  void BoostAlongY(TLorentzVector &vec, double beta_squared) const;

public:
  /// Constructor for a three particle decay with given masses
  RootGenerator(const ComPWA::FourMomentum &CMSP4_,
                const std::vector<double> &FinalStateMasses_, int seed = -1);

  /// Constructor: Information on the decay is obtained from Kinematics
  RootGenerator(
      const Physics::ParticleStateTransitionKinematicsInfo &KinematicsInfo,
      int seed = -1);

  /// Constructor: Information on the decay is provides via lists of initial and
  /// final states
  RootGenerator(std::shared_ptr<PartList> partL, std::vector<pid> finalS,
                std::vector<pid> initialS, int seed = -1);

  virtual ~RootGenerator(){};

  ComPWA::Event generate();

  void setSeed(unsigned int seed);

  unsigned int getSeed() const;

  double uniform(double min, double max);

  double gauss(double mu, double sigma) const;

protected:
  void init();

  TRandom3 UniformRandomGen;

  ComPWA::FourMomentum CMSP4;
  std::vector<double> FinalStateMasses;
  std::vector<TLorentzVector> FinalStateLorentzVectors;
  double MaximumWeight;
  TVector3 CMSBoostVector;
  // total energy in C.M. minus the sum of the masses
  double CMSEnergyMinusMasses;
};

class UniformTwoBodyGenerator : public RootGenerator {
public:
  UniformTwoBodyGenerator(
      const Physics::ParticleStateTransitionKinematicsInfo &KinematicsInfo,
      int seed, double minSq_, double maxSq_)
      : RootGenerator(KinematicsInfo, seed), minSq(minSq_), maxSq(maxSq_) {}
  virtual ComPWA::Event generate();
  virtual UniformTwoBodyGenerator *clone() {
    return (new UniformTwoBodyGenerator(*this));
  }

protected:
  double minSq, maxSq;
};

} // namespace Tools
} // namespace ComPWA

#endif
