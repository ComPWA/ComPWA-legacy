// Copyright (c) 2016, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains PartialAmplitude interface class
///

#ifndef PHYSICS_PartialAmplitude_HPP_
#define PHYSICS_PartialAmplitude_HPP_

#include <vector>
#include <memory>

#include "Core/Parameter.hpp"
#include "Core/Parameter.hpp"
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Physics {

///
/// \class PartialAmplitude
/// PartialAmplitude class is an interface classe that resembles a two-body decay.
///
class PartialAmplitude {

public:
  PartialAmplitude()
      : _preFactor(1, 0), phspVolume_(1), CurrentIntegral(1.0),
        CurrentMagnitude(0.0), CurrentPhase(0.0){};

  virtual ~PartialAmplitude(){};

  virtual PartialAmplitude *clone(std::string newName = "") const = 0;

  //======= INTEGRATION/NORMALIZATION ===========

  /// Get current normalization.
  virtual double normalization() const = 0;

  /// Check of parameters have changed and normalization has to be recalculatecd
  bool isModified() const {
    if (magnitude() != CurrentMagnitude || phase() != CurrentPhase) {
      const_cast<double &>(CurrentMagnitude) = magnitude();
      const_cast<double &>(CurrentPhase) = phase();
      return true;
    }
    return false;
  }

  //================ EVALUATION =================

  /// Value of PartialAmplitude at \param point with normalization factor
  virtual std::complex<double> evaluate(const DataPoint &point) const {
    return evaluateNoNorm(point) * normalization();
  }

  /// Value of PartialAmplitude at \param point without normalization factor
  virtual std::complex<double> evaluateNoNorm(const DataPoint &point) const = 0;

  //============ SET/GET =================

  virtual std::string name() const { return _name; }

  virtual void setName(std::string name) { _name = name; }

  virtual void setPrefactor(std::complex<double> pre) { _preFactor = pre; }

  virtual std::complex<double> prefactor() const { return _preFactor; }

  virtual std::complex<double> coefficient() const {
    return std::polar(magnitude(), phase());
  }

  std::shared_ptr<ComPWA::DoubleParameter> magnitudeParameter() {
    return _magnitude;
  }

  double magnitude() const { return std::fabs(_magnitude->value()); }

  void setMagnitudeParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _magnitude = par;
  }

  void setMagnitude(double par) { _magnitude->setValue(par); }

  std::shared_ptr<ComPWA::DoubleParameter> phaseParameter() {
    return _phase;
  }

  double phase() const { return _phase->value(); }

  void setPhaseParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _phase = par;
  }

  void setPhase(double par) { _phase->setValue(par); }

  virtual void parameters(ParameterList &list) {
    list.AddParameter(magnitudeParameter());
    list.AddParameter(phaseParameter());
  }

  /// Fill vector with parameters (fast). No check is performed if parameters
  /// already exist. \see GetParameters(ParameterList &list)
  virtual void parametersFast(std::vector<double> &list) const {
    list.push_back(magnitude());
    list.push_back(phase());
  }

  /// Update parameters to the values given in \p list
  virtual void updateParameters(const ParameterList &list) = 0;

  /// Set phase space sample.
  /// We use the phase space sample to calculate the normalization. The sample
  /// should be without efficiency applied.
  virtual void
  setPhspSample(std::shared_ptr<std::vector<ComPWA::DataPoint>> phspSample) {
    _phspSample = phspSample;
  };

  virtual void setPhspVolume(double phspVol) { phspVolume_ = phspVol; }

  virtual double phspVolume() const { return phspVolume_; }

  //=========== FUNCTIONTREE =================

  virtual bool hasTree() const { return false; }

  virtual std::shared_ptr<FunctionTree>
  tree(std::shared_ptr<Kinematics> kin, const ComPWA::ParameterList &sample,
          const ComPWA::ParameterList &toySample, std::string suffix) = 0;

protected:
  std::string _name;
  std::shared_ptr<ComPWA::DoubleParameter> _magnitude;
  std::shared_ptr<ComPWA::DoubleParameter> _phase;
  std::complex<double> _preFactor;

  /// Phsp sample for numerical integration
  std::shared_ptr<std::vector<ComPWA::DataPoint>> _phspSample;

  /// The phase-space volume is needed for proper normalization of the PartialAmplitude
  double phspVolume_;

  virtual double integral() const {
    if (!_phspSample->size()) {
      LOG(debug)
          << "CoherentIntensity::Integral() | Integral can not be calculated "
             "since no phsp sample is set. Set a sample using "
             "SetPhspSamples(phspSample, toySample)!";
      return 1.0;
    }

    double sumIntens = 0;
    for (auto i : *_phspSample.get())
      sumIntens += std::norm(evaluateNoNorm(i));

    double integral = (sumIntens * phspVolume_ / _phspSample->size());
    LOG(trace) << "PartialAmplitude::Integral() | Integral is " << integral << ".";
    assert(!std::isnan(integral));
    return integral;
  }

  //! Integral value (temporary)
  double CurrentIntegral;

private:
  //! Temporary value of mass (used to trigger recalculation of normalization)
  double CurrentMagnitude;
  double CurrentPhase;
};

} /* namespace Physics */
} /* namespace ComPWA */
#endif /* PHYSICS_PartialAmplitude_HPP_ */
