
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
  Evaluate(const ComPWA::dataPoint &point) const = 0;

  //============ SET/GET =================
  
  virtual void SetName(std::string n) { _name = n; }

  virtual std::string GetName() { return _name; }

  virtual void SetModified(bool b = true) const {
    const_cast<bool &>(_modified) = b;
    const_cast<double &>(_current_mass) = _mass->GetValue();
  }

  virtual bool GetModified() const { return _modified; }

  virtual void GetParameters(ParameterList &list);


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

  virtual void SetLimits(std::pair<double, double> lim) {
    LOG(trace) << "AbstractDynamicalFunction::SetLimits() | Setting invariant "
                  "mass limits for "
               << _name << " to "
               << "[" << lim.first << "," << lim.second << "].";
    _limits = lim;
  }

  virtual std::pair<double, double> GetLimits() { return _limits; }

  //! Set position of variables within dataPoint
  void SetDataPosition(int pos) { _dataPos = pos; }

  //! Get position of variables within dataPoint
  int GetDataPosition() const { return _dataPos; }

  //! Set position of variables within dataPoint
  void SetSubSystem(SubSystem sys) {
    _dataPos = 3 *
               dynamic_cast<HelicityKinematics *>(Kinematics::Instance())
                   ->GetDataID(sys);
  }

  //=========== FUNCTIONTREE =================
  
  virtual bool HasTree() const { return false; }
  
  /**! Setup function tree */
  virtual std::shared_ptr<ComPWA::FunctionTree>
  GetTree(const ComPWA::ParameterList &sample, std::string suffix = "") = 0;

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

  //! Minimum and Maximum of invariant mass
  std::pair<double, double> _limits;

  //! Resonance spin
  ComPWA::Spin _spin;

  //! Position of invariant mass in dataPoint
  int _dataPos;

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
  std::complex<double> qValue;
  if (spin > 0) {
    qValue = Kinematics::qValue(mR, ma, mb);
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
