// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_RANDOM_H_
#define COMPWA_RANDOM_H_

#include <random>

namespace ComPWA {

///
/// \class UniformRealDistribution
/// \brief Interface class for generating random doubles in the range [0,1)
///
class UniformRealNumberGenerator {
public:
  virtual ~UniformRealNumberGenerator() = default;
  // generate random double in the range from [0,1)
  virtual double operator()() = 0;
  virtual int getSeed() const = 0;
  virtual void setSeed(int seed) = 0;
};

class StdUniformRealGenerator : public UniformRealNumberGenerator {
  std::mt19937 MersenneTwisterRandomGenerator;
  std::uniform_real_distribution<double> UniformDistribution;
  int Seed;

public:
  StdUniformRealGenerator(int seed = 123456);

  double operator()() final;
  int getSeed() const final;
  void setSeed(int seed) final;
};

} // namespace ComPWA

#endif
