/*
 * EvtGenGenerator.hpp
 *
 *  Created on: Nov 21, 2018
 *      Author: steve
 */

#ifndef TOOLS_EVTGENGENERATOR_HPP_
#define TOOLS_EVTGENGENERATOR_HPP_

#include "Core/Generator.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Particle.hpp"
#include "ThirdParty/EvtGen/EvtRandomEngine.hh"

#include <random>

namespace ComPWA {
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
  EvtGenGenerator(std::shared_ptr<PartList> partL,
                  std::shared_ptr<Kinematics> kin, unsigned int seed = -1);

  virtual ~EvtGenGenerator();

  virtual ComPWA::Event generate();

  virtual void setSeed(unsigned int seed);

  virtual unsigned int getSeed() const;

  virtual double uniform(double min, double max);

  virtual double gauss(double mu, double sigma) const { return 0; }
};

} /* namespace Tools */
} /* namespace ComPWA */

#endif /* TOOLS_EVTGENGENERATOR_HPP_ */
