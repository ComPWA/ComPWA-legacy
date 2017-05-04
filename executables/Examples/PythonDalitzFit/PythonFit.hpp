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
/*! \class PythonFit
 * @file PythonFit.hpp
 * This class provides a Manager for simple fits. It creates a set of modules
 * for an unbinned likelihood fit of a three-body final state.
 */

#ifndef _PYTHONFIT_HPP_
#define _PYTHONFIT_HPP_

#include <vector>
#include <memory>

#include "DataReader/Data.hpp"
#include "Estimator/Estimator.hpp"
#include "Core/Amplitude.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/FitResult.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Generator.hpp"

using namespace ComPWA;
using namespace ComPWA::DataReader;

class PythonFit {
public:
  PythonFit();

  virtual ~PythonFit();

  virtual void setData(std::shared_ptr<Data> d) { sampleData_ = d; };
  virtual std::shared_ptr<Data> getData() { return sampleData_; };
  virtual void setBackground(std::shared_ptr<Data> d) { sampleBkg_ = d; };
  virtual std::shared_ptr<Data> getBackground() { return sampleBkg_; };
  virtual void
  setPhspSample(std::shared_ptr<Data> phsp,
                std::shared_ptr<Data> truePhsp = std::shared_ptr<Data>());
  virtual std::shared_ptr<Data> getPhspSample() { return samplePhsp_; };
  virtual void setTruePhspSample(std::shared_ptr<Data>);
  virtual std::shared_ptr<Data> getTruePhspSample() { return sampleTruePhsp_; };

  virtual void setAmplitude(std::shared_ptr<Amplitude> d) { amp_ = d; };
  virtual std::shared_ptr<Amplitude> getAmplitude() { return amp_; };
  virtual void setBkgAmplitude(std::shared_ptr<Amplitude> d) { ampBkg_ = d; };
  virtual std::shared_ptr<Amplitude> getBkgAmplitude() { return ampBkg_; };
  virtual void setOptimizer(std::shared_ptr<Optimizer::Optimizer> d) {
    opti_ = d;
  };
  virtual std::shared_ptr<Optimizer::Optimizer> getOptimizer() {
    return opti_;
  };
  virtual void setGenerator(std::shared_ptr<Generator> d) { gen_ = d; };
  virtual std::shared_ptr<Generator> getGenerator() { return gen_; };

  virtual std::shared_ptr<FitResult> startFit(ParameterList &);

  virtual void StartFit();

  /**Generate phase space events by Hit&Miss
   *
   * @param number Number of events to generate
   * @return
   */
  virtual bool generatePhsp(int number);

  /**Generate signal events by Hit&Miss
   * 1) In case no phsp sample is set and the @param number is larger zero,
   * phsp events are generated on the fly.
   * 2) In case a phsp sample is set and @param number is smaller zero,
   * the whole sample is used for event generation.
   *
   * @param number Number of events to generate
   * @return
   */
  virtual bool generate(int number);

  /**Generate background events by Hit&Miss
   * 1) In case no phsp sample is set and the @param number is larger zero, phsp
   * events
   * are generated on the fly.
   * 2) In case a phsp sample is set and @param number is smaller zero, the
   * whole sample
   * is used for event generation.
   *
   * @param number Number of events to generate
   * @return
   */
  virtual bool generateBkg(int number);

  virtual void SetAmplitudesData(std::vector<std::shared_ptr<Amplitude>> ampVec,
                                 std::vector<double> fraction,
                                 std::vector<std::shared_ptr<Data>> dataVec);

  virtual void GenAmplitudesData(int nEvents);

  virtual std::vector<std::shared_ptr<Data>> GetData() { return _dataVec; }

protected:
  static bool gen(int number, std::shared_ptr<Generator> gen,
                  std::shared_ptr<Amplitude> amp, std::shared_ptr<Data> data,
                  std::shared_ptr<Data> phsp = std::shared_ptr<Data>(),
                  std::shared_ptr<Data> phspTrue = std::shared_ptr<Data>());

  ParameterList par;

  std::shared_ptr<Data> sampleData_; /*!< Pointer to data sample */
  std::shared_ptr<Data> sampleBkg_;  /*!< Pointer to data sample */

  std::shared_ptr<Data> samplePhsp_;     /*!< Pointer to phsp sample */
  std::shared_ptr<Data> sampleTruePhsp_; /*!< Pointer to true phsp sample */

  std::shared_ptr<Amplitude> amp_;    /*!< Pointer to signal model */
  std::shared_ptr<Amplitude> ampBkg_; /*!< Pointer to background model */
  std::shared_ptr<Optimizer::Optimizer>
      opti_;                       /*!< Pointer to Optimizer-Module */
  std::shared_ptr<Generator> gen_; /*!< Pointer to Generator-Module */

  std::vector<std::shared_ptr<Amplitude>> _ampVec;
  std::vector<double> _fraction;
  std::vector<std::shared_ptr<Data>> _dataVec;

  double M; // GeV/c² (J/psi+)
  double Br; // GeV/c² (width)
  double m1; // GeV/c² (gamma)
  double m2; // GeV/c² (pi)
  double m3; // GeV/c² (pi)
  //double c = 299792458.; // m/s
  double PI; // m/s

  unsigned int nFitEvents;
  unsigned int nStartEvent;
};

#endif
