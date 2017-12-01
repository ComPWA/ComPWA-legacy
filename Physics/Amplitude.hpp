// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Physics Interface Base-Class.
/*! \class Amplitude
 * @file Amplitude.hpp
 * This class provides the interface to the model which tries to describe the
 * intensity. As it is pure virtual, one needs at least one implementation to
 * provide an model for the analysis which calculates intensities for an event
 * on
 * basis model parameters. If a new physics-model is derived from and fulfills
 * this base-class, no change in other modules are necessary to work with the
 * new
 * physics module.
 */

#ifndef AMPLITUDE_HPP_
#define AMPLITUDE_HPP_

#include <vector>
#include <memory>
#include <math.h>

#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Physics {

class Amplitude {

public:
  //============ CONSTRUCTION ==================

  //! Constructor with an optional, unique name and an optional efficiency
  Amplitude(std::string name = "")
      : Name(name), PreFactor(1, 0), CurrentMagnitude(0.0),
        CurrentPhase(0.0){};

  //! Destructor
  virtual ~Amplitude() { /* nothing */
  }

  //! Function to create a full copy of the amplitude
  virtual Amplitude *clone(std::string newName = "") const = 0;

  //======= INTEGRATION/NORMALIZATION ===========

  /// Check if parameters have changed and ifnormalization has to be
  /// recalculatecd.
  bool isModified() const {
    if (magnitude() != CurrentMagnitude || phase() != CurrentPhase) {
      const_cast<double &>(CurrentMagnitude) = magnitude();
      const_cast<double &>(CurrentPhase) = phase();
      return true;
    }
    return false;
  }

  //================ EVALUATION =================

  /// Calculate value of amplitude at \p point.
  virtual std::complex<double> evaluate(const DataPoint &point) const = 0;

  //============ SET/GET =================

  virtual std::string name() const { return Name; }

  virtual void setName(std::string name) { Name = name; }

  virtual std::complex<double> coefficient() const {
    return std::polar(magnitude(), phase());
  }

  /// Update parameters to the values given in \p list
  virtual void updateParameters(const ParameterList &par) = 0;

  /// Fill parameters to list
  virtual void parameters(ParameterList &list) {
    list.AddParameter(Magnitude);
    list.AddParameter(Phase);
  }

  /// Fill vector with parameters.
  /// In comparisson to GetParameters(ParameterList &list) no checks are
  /// performed here. So this should be much faster.
  virtual void parametersFast(std::vector<double> &list) const {
    list.push_back(magnitude());
    list.push_back(phase());
  }

  virtual std::shared_ptr<ComPWA::DoubleParameter> magnitudeParameter() {
    return Magnitude;
  }

  virtual double magnitude() const {
    return std::fabs(Magnitude->value());
  }

  virtual void
  setMagnitudeParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    Magnitude = par;
  }

  virtual void setMagnitude(double par) { Magnitude->setValue(par); }

  virtual std::shared_ptr<ComPWA::DoubleParameter> phaseParameter() {
    return Phase;
  }

  virtual double phase() const { return Phase->value(); }

  virtual void setPhaseParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    Phase = par;
  }

  virtual void setPhase(double par) { Phase->setValue(par); }

  virtual void setPreFactor(std::complex<double> pre) { PreFactor = pre; }

  virtual std::complex<double> preFactor() const { return PreFactor; }

  /// Set phase space sample.
  /// We use the phase space sample to calculate the normalization. The sample
  /// should be without efficiency applied.
  virtual void
  setPhspSample(std::shared_ptr<std::vector<ComPWA::DataPoint>> phspSample) = 0;

  //=========== FUNCTIONTREE =================

  virtual bool hasTree() const { return 0; }

  virtual std::shared_ptr<FunctionTree> tree(std::shared_ptr<Kinematics> kin,
                                                const ParameterList &sample,
                                                const ParameterList &toySample,
                                                std::string suffix) = 0;

protected:
  std::string Name;

  std::complex<double> PreFactor;

  std::shared_ptr<DoubleParameter> Magnitude;

  std::shared_ptr<DoubleParameter> Phase;

private:
  double CurrentMagnitude;
  double CurrentPhase;
};

typedef std::vector<std::shared_ptr<Amplitude>>::iterator ampItr;

} // namespace Physics
} // namespace ComPWA
#endif
