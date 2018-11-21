// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef TOOLS_ROOTGENERATOR_HPP_
#define TOOLS_ROOTGENERATOR_HPP_

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "TRandom3.h"

#include "Core/Generator.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Particle.hpp"

namespace ComPWA {
namespace Tools {

class RootGenerator : public Generator {
  double PDK(double a, double b, double c) const;
  void BoostAlongY(TLorentzVector &vec, double beta_squared) const;

public:
  /// Constructor for a three particle decay with given masses
  RootGenerator(const ComPWA::FourMomentum &CMSP4_,
                const std::vector<double> &FinalStateMasses_, int seed = -1);

  /// Constructor: Information on the decay is obtained from Kinematics
  RootGenerator(std::shared_ptr<PartList> partL,
                std::shared_ptr<Kinematics> kin, int seed = -1);

  /// Constructor: Information on the decay is provides via lists of initial and
  /// final states
  RootGenerator(std::shared_ptr<PartList> partL, std::vector<pid> finalS,
                std::vector<pid> initialS, int seed = -1);

  virtual ~RootGenerator(){};

  virtual ComPWA::Event generate();

  virtual void setSeed(unsigned int seed);

  virtual unsigned int getSeed() const;

  virtual double uniform(double min, double max);

  virtual double gauss(double mu, double sigma) const;

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
  UniformTwoBodyGenerator(std::shared_ptr<PartList> partL,
                          std::shared_ptr<Kinematics> kin, int seed,
                          double minSq_, double maxSq_)
      : RootGenerator(partL, kin, seed), minSq(minSq_), maxSq(maxSq_) {
    auto const &KinProps(kin->getKinematicsProperties());
    auto finalS = KinProps.FinalState;
    auto initialS = KinProps.InitialState;
    if (finalS.size() != 2)
      throw std::runtime_error("UniformTwoBodyGenerator::"
                               "UniformTwoBodyGenerator() | Not a two body "
                               "decay!");
    for (auto ParticlePid : finalS) { // particle 0 is mother particle
      FinalStateMasses.push_back(FindParticle(partL, ParticlePid).GetMass());
    }
  }
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
