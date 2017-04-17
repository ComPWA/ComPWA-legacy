 
//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//----------------------------------------------------------------------------------

#ifndef PHYSICS_HELICITYAMPLITUDE_ABSTRACTDYNAMICALFUNCTION_HPP_
#define PHYSICS_HELICITYAMPLITUDE_ABSTRACTDYNAMICALFUNCTION_HPP_

#include <complex>

#include <boost/property_tree/ptree.hpp>

#include "Core/PhysConst.hpp"
#include "Core/DataPoint.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Spin.hpp"
#include "Core/FunctionTree.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

enum normStyle {
  none, /*!< no normaliztion between Amplitudes. */
  /*!< all amplitudes are normalized to one.
   *  The normalization factor is \f$ 1/\sqrt(\int |A|^2)\f$ */
  one
};

class AbstractDynamicalFunction {
public:
  AbstractDynamicalFunction();

  virtual ~AbstractDynamicalFunction();

  virtual std::complex<double> Evaluate(const ComPWA::dataPoint &point,
                                        int pos) const = 0;

  /**! Get current normalization.  */
  virtual double GetNormalization() const = 0;

  void CheckModified() const {
    if (_mass->GetValue() != _current_mass) {
      const_cast<bool &>(_modified) = 1;
      const_cast<double &>(_current_mass) = _mass->GetValue();
    }
    return;
  }

  virtual bool HasTree() const { return false; }

  /**! Setup function tree */
  virtual std::shared_ptr<ComPWA::FunctionTree>
  GetTree(ComPWA::ParameterList &sample,
          ComPWA::ParameterList &toySample, int pos, std::string suffix="") = 0;

  /**
   Set decay width

   @param w Decay width
   */
  virtual void SetMass(std::shared_ptr<DoubleParameter> mass) { _mass = mass; }

  /**
   Get decay mass

   @return Decay mass
   */
  virtual std::shared_ptr<DoubleParameter> GetMass() { return _mass; }

  /**
   Set decay mass

   @param w Decay mass
   */
  virtual void SetMass(double mass) { _mass->SetValue(mass); }

  /**
   Get decay mass

   @return Decay mass
   */
  virtual double GetMassValue() const { return _mass->GetValue(); }

  virtual void SetDecayMassA(double mass) { _massA = mass; }

  virtual double GetDecayMassA() const { return _massA; }

  virtual void SetDecayMassB(double mass) { _massB = mass; }

  virtual double GetDecayMassB() const { return _massB; }

  virtual ComPWA::Spin GetSpin() const { return _spin; }

  virtual void SetSpin(ComPWA::Spin spin) { _spin = spin; }

  virtual void SetModified(bool b = true) const {
    const_cast<bool &>(_modified) = b;
    const_cast<double &>(_current_mass) = _mass->GetValue();
  }

  virtual bool GetModified() const { return _modified; }

protected:
  //! Name of resonance
  std::string _name;

  //! Type of resonance normalization
  normStyle _normStyle;

  //! Precision of MC integration
  int _mcPrecision;

  //! Integral
  virtual double Integral() const { return 1.0; };

  //! Masses of daughter particles
  double _massA, _massB;

  //! Resonance mass
  std::shared_ptr<ComPWA::DoubleParameter> _mass;

  //! Resonance sub system
  unsigned int _dataIndex;

  //! Resonance spin
  ComPWA::Spin _spin;

  //! Integral value (temporary)
  double _current_integral;

private:
  //! Resonance shape was modified (recalculate the normalization)
  bool _modified;

  //! Temporary value of mass (used to trigger recalculation of normalization)
  double _current_mass;
};

inline std::complex<double> widthToCoupling(double mSq, double mR, double width,
                                            double ma, double mb, double spin,
                                            double mesonRadius,
                                            formFactorType type) {
  double sqrtS = sqrt(mSq);

  // calculate gammaA(s_R)
  std::complex<double> gammaA(1, 0); // spin==0
  if (spin > 0) {
    std::complex<double> qValue = Kinematics::qValue(mR, ma, mb);
    double ffR = HelicityKinematics::FormFactor(mR, ma, mb, spin, mesonRadius,
                                                qValue, type);
    std::complex<double> qR = std::pow(qValue, spin);
    gammaA = ffR * qR;
  }

  // calculate phsp factor
  std::complex<double> rho = Kinematics::phspFactor(sqrtS, ma, mb);

  std::complex<double> denom = gammaA * sqrt(rho);
  std::complex<double> res = std::complex<double>(sqrt(mR * width), 0) / denom;

#ifndef NDEBUG
  // check for NaN
  //	if( std::isnan(res.real()) || std::isnan(res.imag()) )
  //		throw
  // std::runtime_error("AmpAbsDynamicalFunction::widthToCoupling()
  //| "
  //				"Result is NaN!");
  // check for inf
  if (std::isinf(res.real()) || std::isinf(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::widthToCoupling() | "
                             "Result is inf!");
#endif

  return res;
}

inline std::complex<double> couplingToWidth(double mSq, double mR, double g,
                                            double ma, double mb, double spin,
                                            double mesonRadius,
                                            formFactorType type,
                                            std::complex<double> phspFactor) {
  // calculate gammaA(s_R)
  std::complex<double> gammaA(1, 0); // spin==0

  if (spin > 0 || type == formFactorType::CrystalBarrel) {
    std::complex<double> qValue = Kinematics::qValue(mR, ma, mb);
    double ffR = HelicityKinematics::FormFactor(mR, ma, mb, spin, mesonRadius,
                                                qValue, type);
    std::complex<double> qR = std::pow(qValue, spin);
    gammaA = ffR * qR;
  }

  // calculate phsp factor
  std::complex<double> res = std::norm(gammaA) * g * g * phspFactor / mR;

#ifndef NDEBUG
  // check for NaN
  if (std::isnan(res.real()) || std::isnan(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::couplingToWidth() | "
                             "Result is NaN!");
  // check for inf
  if (std::isinf(res.real()) || std::isinf(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::couplingToWidth() | "
                             "Result is inf!");
#endif

  return res;
}

inline std::complex<double> couplingToWidth(double mSq, double mR, double g,
                                            double ma, double mb, double spin,
                                            double mesonRadius,
                                            formFactorType type) {
  double sqrtM = sqrt(mSq);
  std::complex<double> phspFactor = Kinematics::phspFactor(sqrtM, ma, mb);

  return couplingToWidth(mSq, mR, g, ma, mb, spin, mesonRadius, type,
                         phspFactor);
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_ABSTRACTDYNAMICALFUNCTION_HPP_ */
