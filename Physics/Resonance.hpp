/*
 * Resonance.hpp
 *
 *  Created on: Mar 3, 2016
 *      Author: weidenka
 */

#ifndef CORE_RESONANCE_HPP_
#define CORE_RESONANCE_HPP_

#include <vector>
#include <memory>

#include <boost/iterator/filter_iterator.hpp>

#include "Core/AbsParameter.hpp"
#include "Core/Parameter.hpp"
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"

namespace ComPWA {
namespace Physics {

enum normStyle {
  none, /*!< no normaliztion between Amplitudes. */
  /*!< all amplitudes are normalized to one.
   *  The normalization factor is \f$ 1/\sqrt(\int |A|^2)\f$ */
  one
};

class Resonance {

public:
  //============ CONSTRUCTION ==================

  Resonance() : _preFactor(1, 0){};

  virtual ~Resonance(){};

  //! Clone function
  virtual Resonance *Clone(std::string newName = "") const = 0;

  //======= INTEGRATION/NORMALIZATION ===========

  /**! Get current normalization.  */
  virtual double GetNormalization() const = 0;

  //! Check of parameters have changed and normalization has to be recalculatecd
  bool CheckModified() const {
    if (GetMagnitude() != _current_magnitude || GetPhase() != _current_phase) {
      const_cast<double &>(_current_magnitude) = GetMagnitude();
      const_cast<double &>(_current_phase) = GetPhase();
      return true;
    }
    return false;
  }

  //================ EVALUATION =================

  //! value of resonance at \param point
  virtual std::complex<double> Evaluate(const dataPoint &point) const {
    return EvaluateNoNorm(point) * GetNormalization();
  }

  //! value of resonance at \param point
  virtual std::complex<double> EvaluateNoNorm(const dataPoint &point) const = 0;

  //============ SET/GET =================

  //! Get resonance name
  virtual std::string GetName() const { return _name; }

  //! Set resonance name
  virtual void SetName(std::string name) { _name = name; }

  //! Set prefactor
  virtual void SetPrefactor(std::complex<double> pre) { _preFactor = pre; }

  //! Get prefactor
  virtual std::complex<double> GetPrefactor() const { return _preFactor; }

  //! Get coefficient
  virtual std::complex<double> GetCoefficient() const {
    return std::polar(GetMagnitude(), GetPhase());
  }

  /**
   Get strength parameter

   @return strength parameter
   */
  std::shared_ptr<ComPWA::DoubleParameter> GetMagnitudeParameter() {
    return _magnitude;
  }

  /**
   Get strength parameter

   @return strength parameter
   */
  double GetMagnitude() const { return std::fabs(_magnitude->GetValue()); }

  /**
   Set strength parameter

   @param par Strength parameter
   */
  void SetMagnitudeParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _magnitude = par;
  }

  /**
   Set strength parameter

   @param par Strength parameter
   */
  void SetMagnitude(double par) { _magnitude->SetValue(par); }

  /**
   Get phase parameter

   @return Phase parameter
   */
  std::shared_ptr<ComPWA::DoubleParameter> GetPhaseParameter() {
    return _phase;
  }

  /**
   Get phase parameter

   @return Phase parameter
   */
  double GetPhase() const { return _phase->GetValue(); }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  void SetPhaseParameter(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _phase = par;
  }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  void SetPhase(double par) { _phase->SetValue(par); }

  virtual void GetParameters(ParameterList &list) {
    list.AddParameter(GetMagnitudeParameter());
    list.AddParameter(GetPhaseParameter());
  }

  /*! Set phase space sample
   * We use the phase space sample to calculate the normalization. The sample
   * should be without efficiency applied.
   */
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample) {
    _phspSample = phspSample;
  };

  //=========== FUNCTIONTREE =================

  //! Check of tree is available
  virtual bool HasTree() const { return false; }

  virtual std::shared_ptr<FunctionTree>
  GetTree(const ComPWA::ParameterList &sample,
          const ComPWA::ParameterList &toySample, std::string suffix) = 0;

protected:
  std::string _name;
  std::shared_ptr<ComPWA::DoubleParameter> _magnitude;
  std::shared_ptr<ComPWA::DoubleParameter> _phase;
  std::complex<double> _preFactor;

  //! Phsp sample for numerical integration
  std::shared_ptr<std::vector<ComPWA::dataPoint>> _phspSample;

  //! Integral
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

    double phspVol = Kinematics::Instance()->GetPhspVolume();
    double integral = (sumIntens * phspVol / _phspSample->size());
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
#endif /* CORE_RESONANCE_HPP_ */
