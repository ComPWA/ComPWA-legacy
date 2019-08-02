// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Random.hpp"

namespace ComPWA {

StdUniformRealGenerator::StdUniformRealGenerator(int seed)
    : MersenneTwisterRandomGenerator(seed), UniformDistribution(0.0, 1.0),
      Seed(seed) {}

double StdUniformRealGenerator::operator()() {
  return UniformDistribution(MersenneTwisterRandomGenerator);
}

int StdUniformRealGenerator::getSeed() const { return Seed; }

void StdUniformRealGenerator::setSeed(int seed) {
  Seed = seed;
  MersenneTwisterRandomGenerator.seed(seed);
}

} /* namespace ComPWA */
