// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef TOOLS_ROOTGENERATOR_HPP_
#define TOOLS_ROOTGENERATOR_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

#include "TGenPhaseSpace.h"

#include "Core/Generator.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Tools {

class RootGenerator : public Generator {

public:
  /// Constructor for a three particle decay with given masses
  RootGenerator(double sqrtS, double m1, double m2, double m3, int seed = -1);

  /// Constructor: Information on the decay is obtained from Kinematics
  RootGenerator(std::shared_ptr<PartList> partL,
                std::shared_ptr<Kinematics> kin, int seed = -1);

  /// Constructor: Information on the decay is provides via lists of initial and
  /// final states
  RootGenerator(std::shared_ptr<PartList> partL, std::vector<pid> finalS,
                std::vector<pid> initialS, int seed = -1);

  ~RootGenerator() { delete[] masses; };

  virtual RootGenerator *clone();

  virtual void generate(Event &evt);

  virtual void setSeed(unsigned int seed);

  virtual unsigned int seed() const;

  virtual double uniform(double min, double max) const;

  virtual double gauss(double mu, double sigma) const;

  virtual TGenPhaseSpace *GetGenerator() { return &event; }

protected:
  TGenPhaseSpace event;

  size_t nPart;

  Double_t *masses;
  FourMomentum cmsP4;
};

class UniformTwoBodyGenerator : public RootGenerator {
public:
  UniformTwoBodyGenerator(std::shared_ptr<PartList> partL,
                          std::shared_ptr<Kinematics> kin, int seed,
                          double minSq_, double maxSq_)
      : RootGenerator(partL, kin, seed), minSq(minSq_), maxSq(maxSq_) {
    if (kin->getKinematicsProperties().FinalState.size() != 2)
      throw std::runtime_error("UniformTwoBodyGenerator::"
                               "UniformTwoBodyGenerator() | Not a two body "
                               "decay!");
  }
  virtual void generate(Event &evt);
  virtual UniformTwoBodyGenerator *clone() {
    return (new UniformTwoBodyGenerator(*this));
  }

protected:
  double minSq, maxSq;
};

} // ns::Tools
} // ns::ComPWA

#endif
