//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding functionality to generate set of
//events
//-------------------------------------------------------------------------------
//! Run-Manager for a simple fit.
/*! \class RunManager
 * @file RunManager.hpp
 * This class provides a RunManager for simple fits. To use it, you create
 * all modules you want to use and provide them to the RunManger. It checks
 * for compatibility and if set up correctly it starts the fitting procedure.
 */

/*! \mainpage ComPWA Documentation
 *
 * \section about About ComPWA
 *
 * ComPWA aims to provide a framework and toolkit for amplitude analysis in
 * a very general approach with clearly defined interfaces to extend its
 * features. Therefore, it can easily be extended for different experiments
 * and various models while maintaining a common structure, producing
 * comparable output, using the same fit-procedures and quality assurance
 * tools.
 *
 * \subsection structure Overview
 *
 * The basic modularization of the software can be seen in the following
 * picture. The four basic modules (Experiment, Physics, Estimation and
 * Optimization) each provide interfaces which need to be implemented
 * depending on your problem. They are represented in the framework by the
 * folder structure as well as in the build libraries. When attacking a new
 * problem using ComPWA, you might find that you can use existing module
 * implementations and only need to rewrite or extend part of it.
 *
 * \image html ComPWA_Modules.jpg
 * \image latex ComPWA_Modules.eps "ComPWA Modules" width=10cm
 *
 * \section install Installation
 *
 * \subsection step1 Build Tool: Boost.Build
 *
 * The ComPWA libraries and example-executables are build using the build tool
 * provided by the boost library. First, one needs to provide some paths to
 * external tools via environmental libraries, as can be seen in the setEnv
 * blank file. Secondly, you call the configure.pl which looks for available
 * modules and asks which you want to build. Now you can call bjam to compile
 * the code. All libraries will be located in the lib folder, all binaries in
 * the bin folder. Boost.Build uses Jamfiles for the build configuration. The
 * main file used when bjam is called is Jamroot in the root directory of
 * ComPWA. In the subdirectories of the interfaces and the actual
 * implementations you can find Jamfiles to manage the build process.
 *
 * \subsection step2 Used External Libraries
 *
 * The core Framework only needs a boost installation and a compiler
 * supporting c++11. But for the examples and most likely also for starting
 * some fits and plotting some output, it is recommended to have a root
 * installation (for plots and maybe DataReader::Data storage/reading) and a
 * Minuit2
 * (for optimization) installation ready. Besides Minuit2 also Geneva can be
 * used as optimizer if a compatible installation is available.
 *
 * \section starting Getting Started
 * A good point to start learning how ComPWA works is by looking in the
 * IntegrationTestApp.cpp executable, where a simple fit is performed using
 * all modules of the framework. ...
 */

#ifndef _RUNMANAGER_HPP_
#define _RUNMANAGER_HPP_

#include <vector>
#include <memory>

#include "DataReader/Data.hpp"
#include "Estimator/Estimator.hpp"
#include "Core/AmpIntensity.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/FitResult.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Generator.hpp"

#include "Tools/FitFractions.hpp"
namespace ComPWA {

using namespace DataReader;

class RunManager {
  
public:
  RunManager();
  
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

  virtual void SetAmplitude(std::shared_ptr<AmpIntensity> intens) { intens_ = intens; };
  
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
  virtual bool Generate(int number);

protected:
  static bool gen(int number, std::shared_ptr<Generator> gen,
                  std::shared_ptr<AmpIntensity> amp, std::shared_ptr<Data> data,
                  std::shared_ptr<Data> phsp = std::shared_ptr<Data>(),
                  std::shared_ptr<Data> phspTrue = std::shared_ptr<Data>());

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
