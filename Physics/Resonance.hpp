// Copyright (c) 2016, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains Resonance interface class
///

#ifndef PHYSICS_RESONANCE_HPP_
#define PHYSICS_RESONANCE_HPP_

#include <vector>
#include <memory>

#include "Core/AbsParameter.hpp"
#include "Core/Parameter.hpp"
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Physics {

///
/// \class Resonance
/// Resonance class is an interface classe that resembles a two-body decay.
///
class Resonance {

public:
  Resonance()
      : _preFactor(1, 0), phspVolume_(1), _current_integral(1.0),
        _current_magnitude(0.0), _current_phase(0.0){};

  virtual ~Resonance(){};

  virtual Resonance *Clone(std::string newName = "") const = 0;

  //======= INTEGRATION/NORMALIZATION ===========

  /// Get current normalization.
  virtual double GetNormalization() const = 0;

  /// Check of parameters have changed and normalization has to be recalculatecd
  bool CheckModified() const {
    if (GetMagnitude() != _current_magnitude || GetPhase() != _current_phase) {
      const_cast<double &>(_current_magnitude) = GetMagnitude();
      const_cast<double &>(_current_phase) = GetPhase();
      return true;
    }
    return false;
  }

  //================ EVALUATION =================

  /// Value of resonance at \param point with normalization factor
  virtual std::complex<double> Evaluate(const dataPoint &point) const {
    return EvaluateNoNorm(point) * GetNormalization();
  }

  /// Value of resonance at \param point without normalization factor
  virtual std::complex<double> EvaluateNoNorm(const dataPoint &point) const = 0;

  //============ SET/GET =================

  virtual std::string GetName() const { return _name; }

  virtual void SetName(std::string name) { _name = name; }

  virtual void SetPrefactor(std::complex<double> pre) { _preFactor = pre; }

  virtual std::complex<double> GetPrefactor() const { return _preFactor; }

  virtual std::complex<double> GetCoefficient() const {
    return std::polar(GetMagnitude(), GetPhase());
  }

  std::shared_ptr<ComPWA::DoubleParameter> GetMagnitudeParameter() {
    return _magnitude;
  }

  double GetMagnitude() const { return std::fabs(_magnitude->GetValue()); }

  void SetMagnitudeParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _magnitude = par;
  }

  void SetMagnitude(double par) { _magnitude->SetValue(par); }

  std::shared_ptr<ComPWA::DoubleParameter> GetPhaseParameter() {
    return _phase;
  }

  double GetPhase() const { return _phase->GetValue(); }

  void SetPhaseParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _phase = par;
  }

  void SetPhase(double par) { _phase->SetValue(par); }

  virtual void GetParameters(ParameterList &list) {
    list.AddParameter(GetMagnitudeParameter());
    list.AddParameter(GetPhaseParameter());
  }

  /// Fill vector with parameters (fast). No check is performed if parameters
  /// already exist. \see GetParameters(ParameterList &list)
  virtual void GetParametersFast(std::vector<double> &list) const {
    list.push_back(GetMagnitude());
    list.push_back(GetPhase());
  }

  /// Update parameters to the values given in \p list
  virtual void UpdateParameters(const ParameterList &list) = 0;

  /// Set phase space sample.
  /// We use the phase space sample to calculate the normalization. The sample
  /// should be without efficiency applied.
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample) {
    _phspSample = phspSample;
  };

  virtual void SetPhspVolume(double phspVol) { phspVolume_ = phspVol; }

  virtual double GetPhspVolume() const { return phspVolume_; }

  //=========== FUNCTIONTREE =================

  virtual bool HasTree() const { return false; }

  virtual std::shared_ptr<FunctionTree>
  GetTree(std::shared_ptr<Kinematics> kin, const ComPWA::ParameterList &sample,
          const ComPWA::ParameterList &toySample, std::string suffix) = 0;

protected:
  std::string _name;
  std::shared_ptr<ComPWA::DoubleParameter> _magnitude;
  std::shared_ptr<ComPWA::DoubleParameter> _phase;
  std::complex<double> _preFactor;

  /// Phsp sample for numerical integration
  std::shared_ptr<std::vector<ComPWA::dataPoint>> _phspSample;

  /// The phase-space volume is needed for proper normalization of the resonance
  double phspVolume_;

  virtual double Integral() const {
    if (!_phspSample->size()) {
      LOG(debug)
          << "CoherentIntensity::Integral() | Integral can not be calculated "
             "since no phsp sample is set. Set a sample using "
             "SetPhspSamples(phspSample, toySample)!";
      return 1.0;
    }

    double sumIntens = 0;
    for (auto i : *_phspSample.get())
      sumIntens += std::norm(EvaluateNoNorm(i));

    double integral = (sumIntens * phspVolume_ / _phspSample->size());
    LOG(trace) << "Resonance::Integral() | Integral is " << integral << ".";
    assert(!std::isnan(integral));
    return integral;
  }

  //! Integral value (temporary)
  double _current_integral;

private:
  //! Temporary value of mass (used to trigger recalculation of normalization)
  double _current_magnitude;
  double _current_phase;
};

} /* namespace Physics */
} /* namespace ComPWA */
#endif /* PHYSICS_RESONANCE_HPP_ */
