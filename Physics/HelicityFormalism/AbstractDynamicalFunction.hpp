
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

class AbstractDynamicalFunction {
  
public:
  //============ CONSTRUCTION ==================
  
  AbstractDynamicalFunction()
      : _daughterMasses(std::pair<double, double>(-999, -999)),
        _current_mass(-999){};

  virtual ~AbstractDynamicalFunction(){};

  //======= INTEGRATION/NORMALIZATION ===========
  
  bool CheckModified() const {
    if (GetMass() != _current_mass) {
      SetModified();
      const_cast<double &>(_current_mass) = _mass->GetValue();
      return true;
    }
    return false;
  }

  //================ EVALUATION =================
  
  virtual std::complex<double>
  Evaluate(const ComPWA::dataPoint &point, int pos) const = 0;

  //============ SET/GET =================
  
  virtual void SetName(std::string n) { _name = n; }

  virtual std::string GetName() { return _name; }

  virtual void SetModified(bool b = true) const {
    const_cast<bool &>(_modified) = b;
    const_cast<double &>(_current_mass) = _mass->GetValue();
  }

  virtual bool GetModified() const { return _modified; }

  virtual void GetParameters(ParameterList &list);

  //! Fill vector with parameters
  virtual void GetParametersFast(std::vector<double> &list) const {
    list.push_back(GetMass());
  }

  /**
   Set decay mass

   @param w Decay mass
   */
  virtual void SetMassParameter(std::shared_ptr<DoubleParameter> mass) {
    _mass = mass;
  }

  /**
   Get decay mass

   @return Decay mass
   */
  virtual std::shared_ptr<DoubleParameter> GetMassParameter() { return _mass; }

  /**
   Set decay mass

   @param w Decay mass
   */
  virtual void SetMass(double mass) { _mass->SetValue(mass); }

  /**
   Get decay mass

   @return Decay mass
   */
  virtual double GetMass() const { return _mass->GetValue(); }

  virtual void SetDecayMasses(std::pair<double, double> m) {
    _daughterMasses = m;
  }

  virtual std::pair<double, double> GetDecayMasses() const {
    return _daughterMasses;
  }

  virtual void SetDecayNames(std::pair<std::string, std::string> n) {
    _daughterNames = n;
  }

  virtual std::pair<std::string, std::string> GetDecayNames() const {
    return _daughterNames;
  }

  virtual ComPWA::Spin GetSpin() const { return _spin; }

  virtual void SetSpin(ComPWA::Spin spin) { _spin = spin; }

  //=========== FUNCTIONTREE =================
  
  virtual bool HasTree() const { return false; }
  
  /**! Setup function tree */
  virtual std::shared_ptr<ComPWA::FunctionTree>
  GetTree(const ComPWA::ParameterList &sample, int pos, std::string suffix = "") = 0;

protected:
  //! Name of resonance
  std::string _name;

  //! Precision of MC integration
  int _mcPrecision;

  //! Masses of daughter particles
  std::pair<double, double> _daughterMasses;

  //! Names of daughter particles
  std::pair<std::string, std::string> _daughterNames;

  //! Resonance mass
  std::shared_ptr<ComPWA::DoubleParameter> _mass;

  //! Resonance spin
  ComPWA::Spin _spin;

private:
  //! Resonance shape was modified (recalculate the normalization)
  bool _modified;

  //! Temporary value of mass (used to trigger recalculation of normalization)
  double _current_mass;
};

/// Calculate form factor
inline double FormFactor(double sqrtS, double ma, double mb,
                                      double spin, double mesonRadius,
                                      std::complex<double> qValue,
                                      formFactorType type) {
  if (mesonRadius == 0)
    return 1.0; // disable form factors
  if (type == formFactorType::noFormFactor)
    return 1.0; // disable form factors
  if (type == formFactorType::BlattWeisskopf && spin == 0) {
    return 1.0;
  }

  // Form factor for a0(980) used by Crystal Barrel (Phys.Rev.D78-074023)
  if (type == formFactorType::CrystalBarrel) {
    if (spin == 0) {
      double qSq = std::norm(qValue);
      double alpha = mesonRadius * mesonRadius / 6;
      return std::exp(-alpha * qSq);
    } else
      throw std::runtime_error("Kinematics::FormFactor() | "
                               "Form factors of type " +
                               std::string(formFactorTypeString[type]) +
                               " are implemented for spin 0 only!");
  }

  // Blatt-Weisskopt form factors with normalization F(x=mR) = 1.
  // Reference: S.U.Chung Annalen der Physik 4(1995) 404-430
  // z = q / (interaction range). For the interaction range we assume
  // 1/mesonRadius
  if (type == formFactorType::BlattWeisskopf) {
    if (spin == 0)
      return 1;
    double qSq = std::norm(qValue);
    double z = qSq * mesonRadius * mesonRadius;
    /* Events below threshold
     * What should we do if event is below threshold? Shouldn't really influence
     * the result
     * because resonances at threshold don't have spin(?) */
    z = std::fabs(z);

    if (spin == 1) {
      return (sqrt(2 * z / (z + 1)));
    } else if (spin == 2) {
      return (sqrt(13 * z * z / ((z - 3) * (z - 3) + 9 * z)));
    } else if (spin == 3) {
      return (
          sqrt(277 * z * z * z / (z * (z - 15) * (z - 15) + 9 * (2 * z - 5))));
    } else if (spin == 4) {
      return (sqrt(12746 * z * z * z * z /
                   ((z * z - 45 * z + 105) * (z * z - 45 * z + 105) +
                    25 * z * (2 * z - 21) * (2 * z - 21))));
    } else
      throw std::runtime_error(
          "Kinematics::FormFactor() | Form factors of type " +
          std::string(formFactorTypeString[type]) +
          " are implemented for spins up to 4!");
  }
  throw std::runtime_error("Kinematics::FormFactor() | Form factor type " +
                           std::to_string((long long int)type) +
                           " not specified!");

  return 0;
}


/// Calculate form factor
inline double FormFactor(double sqrtS, double ma, double mb,
                                      double spin, double mesonRadius,
                                      formFactorType type) {
  if (type == formFactorType::noFormFactor) {
    return 1.0;
  }
  if (type == formFactorType::BlattWeisskopf && spin == 0) {
    return 1.0;
  }

  std::complex<double> qValue = Kinematics::qValue(sqrtS, ma, mb);
  return FormFactor(sqrtS, ma, mb, spin, mesonRadius,
                                        qValue, type);
}

inline std::complex<double> widthToCoupling(double mSq, double mR, double width,
                                            double ma, double mb, double spin,
                                            double mesonRadius,
                                            formFactorType type) {
  double sqrtS = sqrt(mSq);

  // calculate gammaA(s_R)
  std::complex<double> gammaA(1, 0); // spin==0
  std::complex<double> qValue;
  if (spin > 0) {
    qValue = Kinematics::qValue(mR, ma, mb);
    double ffR = FormFactor(mR, ma, mb, spin, mesonRadius,
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
    throw std::runtime_error("AbstractDynamicalFunction::widthToCoupling() | "
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
    double ffR = FormFactor(mR, ma, mb, spin, mesonRadius,
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
