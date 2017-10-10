// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains the RunManger class which controls a simple analysis with ComPWA.
///

#ifndef _RUNMANAGER_HPP_
#define _RUNMANAGER_HPP_

#include <vector>
#include <memory>

#include "Core/Estimator.hpp"
#include "Core/AmpIntensity.hpp"
#include "Core/FitResult.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Generator.hpp"
#include "DataReader/Data.hpp"
#include "Optimizer/Optimizer.hpp"

#include "Tools/FitFractions.hpp"

namespace ComPWA {

using namespace DataReader;

/// \class RunManager
/// Run-Manager for a simple fit.
/// This class provides a RunManager for simple fits. To use it, you create
/// all modules you want to use and provide them to the RunManger. It checks
/// for compatibility and if set up correctly it starts the fitting procedure.
class RunManager {

public:
  RunManager(){};

  RunManager(std::shared_ptr<DataReader::Data>, std::shared_ptr<AmpIntensity>,
             std::shared_ptr<Optimizer::Optimizer>); // Fit

  RunManager(unsigned int size, std::shared_ptr<AmpIntensity>,
             std::shared_ptr<Generator>); // Generate

  virtual ~RunManager();

  virtual void SetData(std::shared_ptr<Data> d) { sampleData_ = d; };

  virtual std::shared_ptr<Data> GetData() { return sampleData_; };

  virtual void
  SetPhspSample(std::shared_ptr<Data> phsp,
                std::shared_ptr<Data> truePhsp = std::shared_ptr<Data>());

  virtual std::shared_ptr<Data> GetPhspSample() { return samplePhsp_; };

  virtual void SetTruePhspSample(std::shared_ptr<Data>);

  virtual std::shared_ptr<Data> GetTruePhspSample() { return sampleTruePhsp_; };

  virtual void SetAmplitude(std::shared_ptr<AmpIntensity> intens) {
    intens_ = intens;
  };

  virtual std::shared_ptr<AmpIntensity> GetAmplitude() { return intens_; };

  virtual void SetOptimizer(std::shared_ptr<Optimizer::Optimizer> d) {
    opti_ = d;
  };

  virtual std::shared_ptr<Optimizer::Optimizer> GetOptimizer() {
    return opti_;
  };

  virtual void SetGenerator(std::shared_ptr<Generator> d) { gen_ = d; };

  virtual std::shared_ptr<Generator> GetGenerator() { return gen_; };

  virtual std::shared_ptr<FitResult> Fit(ParameterList &);

  /**Generate phase space events by Hit&Miss
   *
   * @param number Number of events to generate
   * @return
   */
  virtual bool GeneratePhsp(int number);

  /**Generate signal events by Hit&Miss
   * 1) In case no phsp sample is set and the @param number is larger zero,
   * phsp events are generated on the fly.
   * 2) In case a phsp sample is set and @param number is smaller zero,
   * the whole sample is used for event generation.
   *
   * @param number Number of events to generate
   * @return
   */
  virtual bool Generate(std::shared_ptr<Kinematics> kin, int nEvents);

  static bool gen(int number, std::shared_ptr<Kinematics> kin,
                  std::shared_ptr<Generator> gen,
                  std::shared_ptr<AmpIntensity> amp, std::shared_ptr<Data> data,
                  std::shared_ptr<Data> phsp = std::shared_ptr<Data>(),
                  std::shared_ptr<Data> phspTrue = std::shared_ptr<Data>());

  static bool genPhsp(int nEvents, std::shared_ptr<Generator> gen,
                      std::shared_ptr<Data> sample);

protected:
  std::shared_ptr<Data> sampleData_; /*!< Pointer to data sample */

  std::shared_ptr<Data> samplePhsp_;     /*!< Pointer to phsp sample */
  std::shared_ptr<Data> sampleTruePhsp_; /*!< Pointer to true phsp sample */

  std::shared_ptr<Optimizer::Optimizer>
      opti_;                       /*!< Pointer to Optimizer-Module */
  std::shared_ptr<Generator> gen_; /*!< Pointer to Generator-Module */

  std::shared_ptr<AmpIntensity> intens_;
};

} /* namespace ComPWA */

#endif
