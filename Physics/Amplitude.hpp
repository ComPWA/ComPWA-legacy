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
      : _name(name), _preFactor(1, 0), _current_magnitude(0.0),
        _current_phase(0.0){};

  //! Destructor
  virtual ~Amplitude() { /* nothing */
  }

  //! Function to create a full copy of the amplitude
  virtual Amplitude *Clone(std::string newName = "") const = 0;

  //======= INTEGRATION/NORMALIZATION ===========

  /// Check if parameters have changed and ifnormalization has to be
  /// recalculatecd.
  bool CheckModified() const {
    if (GetMagnitude() != _current_magnitude || GetPhase() != _current_phase) {
      const_cast<double &>(_current_magnitude) = GetMagnitude();
      const_cast<double &>(_current_phase) = GetPhase();
      return true;
    }
    return false;
  }

  //================ EVALUATION =================

  /// Calculate value of amplitude at \p point.
  virtual std::complex<double> Evaluate(const dataPoint &point) const = 0;

  //============ SET/GET =================

  virtual std::string GetName() const { return _name; }

  virtual void SetName(std::string name) { _name = name; }

  virtual std::complex<double> GetCoefficient() const {
    return std::polar(GetMagnitude(), GetPhase());
  }

  /// Update parameters to the values given in \p list
  virtual void UpdateParameters(const ParameterList &par) = 0;

  /// Fill parameters to list
  virtual void GetParameters(ParameterList &list) {
    std::shared_ptr<DoubleParameter> tmp;
    try { // catch BadParameter
      tmp = list.GetDoubleParameter(_phase->GetName());
      // catch and throw std::runtime_error due to failed parameter comparisson
      try {
        if (*tmp == *_phase)
          _phase = tmp;
      } catch (std::exception &ex) {
        throw;
      }
    } catch (BadParameter &ex) {
      list.AddParameter(_phase);
    }

    try { // catch BadParameter
      tmp = list.GetDoubleParameter(_magnitude->GetName());
      // catch and throw std::runtime_error due to failed parameter comparisson
      try {
        if (*tmp == *_magnitude)
          _magnitude = tmp;
      } catch (std::exception &ex) {
        throw;
      }
    } catch (BadParameter &ex) {
      list.AddParameter(_magnitude);
    }
  }

  /// Fill vector with parameters.
  /// In comparisson to GetParameters(ParameterList &list) no checks are
  /// performed here. So this should be much faster.
  virtual void GetParametersFast(std::vector<double> &list) const {
    list.push_back(GetMagnitude());
    list.push_back(GetPhase());
  }

  virtual std::shared_ptr<ComPWA::DoubleParameter> GetMagnitudeParameter() {
    return _magnitude;
  }

  virtual double GetMagnitude() const {
    return std::fabs(_magnitude->GetValue());
  }

  virtual void
  SetMagnitudeParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _magnitude = par;
  }

  virtual void SetMagnitude(double par) { _magnitude->SetValue(par); }

  virtual std::shared_ptr<ComPWA::DoubleParameter> GetPhaseParameter() {
    return _phase;
  }

  virtual double GetPhase() const { return _phase->GetValue(); }

  virtual void SetPhaseParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _phase = par;
  }

  virtual void SetPhase(double par) { _phase->SetValue(par); }

  virtual void SetPreFactor(std::complex<double> pre) { _preFactor = pre; }

  virtual std::complex<double> GetPreFactor() const { return _preFactor; }

  /// Set phase space sample.
  /// We use the phase space sample to calculate the normalization. The sample
  /// should be without efficiency applied.
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample) = 0;

  //=========== FUNCTIONTREE =================

  virtual bool HasTree() const { return 0; }

  virtual std::shared_ptr<FunctionTree> GetTree(std::shared_ptr<Kinematics> kin,
                                                const ParameterList &sample,
                                                const ParameterList &toySample,
                                                std::string suffix) = 0;

protected:
  std::string _name;

  std::complex<double> _preFactor;

  std::shared_ptr<DoubleParameter> _magnitude;

  std::shared_ptr<DoubleParameter> _phase;

private:
  double _current_magnitude;
  double _current_phase;
};

typedef std::vector<std::shared_ptr<Amplitude>>::iterator ampItr;

} /* namespace Physics */
} /* namespace ComPWA */
#endif
