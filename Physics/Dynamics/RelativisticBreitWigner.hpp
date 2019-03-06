// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_DYNAMICS_RELATIVISTICBREITWIGNER_HPP_
#define PHYSICS_DYNAMICS_RELATIVISTICBREITWIGNER_HPP_

#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <vector>

#include "AbstractDynamicalFunction.hpp"
#include "Core/Exceptions.hpp"
#include "Core/Functions.hpp"
#include "Core/Spin.hpp"
#include "FormFactor.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

/// \class RelativisticBreitWigner
/// Relativistic Breit-Wigner model with barrier factors.
/// The dynamical function implemented here is taken from PDG2014 (Eq.47.22)
/// for the one channel case. The dynamic reaction
/// \f[
/// \mathcal{A}_R(s) = \frac{g_p*g}{s - M_R^2 + i \Gamma_R B^2}
/// \f]
/// \f$ g_p, g\f$ are the coupling constants for production and decay and
/// the barrier term \f$ B^2\f$ is parametrised according to Eq.47.23:
/// \f[
///     B^2 = \left( \frac{q(s)}{q(s_R)} \right)^{2L+1} \times
///                     \left( \frac{F(s)}{F(s_R)} \right)^{2}
/// \f]
///
/// NOTE:
/// Make sure to update the orbital angular momentum after construction,
/// otherwise spin J will be used (which is fine is simple cases).
class RelativisticBreitWigner
    : public ComPWA::Physics::Dynamics::AbstractDynamicalFunction {

public:
  //============ CONSTRUCTION ==================
  RelativisticBreitWigner(std::string name,
                          std::pair<std::string, std::string> daughters,
                          std::shared_ptr<ComPWA::PartList> partL);
  virtual ~RelativisticBreitWigner();

  boost::property_tree::ptree save() const;

  //================ EVALUATION =================

  std::complex<double> evaluate(const ComPWA::DataPoint &point,
                                unsigned int pos) const;

  /// Dynamical Breit-Wigner function.
  /// \param mSq Invariant mass squared
  /// \param mR Mass of the resonant state
  /// \param ma Mass of daughter particle
  /// \param mb Mass of daughter particle
  /// \param width Decay width
  /// \param L Orbital angular momentum between two daughters a and b
  /// \param mesonRadius Meson Radius
  /// \param ffType Form factor type
  /// \return Amplitude value
  static std::complex<double>
  dynamicalFunction(double mSq, double mR, double ma, double mb, double width,
                    unsigned int L, double mesonRadius, FormFactorType ffType);

  //============ SET/GET =================

  void SetWidthParameter(std::shared_ptr<ComPWA::FitParameter> w) { Width = w; }

  std::shared_ptr<ComPWA::FitParameter> GetWidthParameter() { return Width; }

  void SetWidth(double w) { Width->setValue(w); }

  double GetWidth() const { return Width->value(); }

  void SetOrbitalAngularMomentum(const ComPWA::Spin &L_) { L = L_; }

  void SetMesonRadiusParameter(std::shared_ptr<ComPWA::FitParameter> r) {
    MesonRadius = r;
  }

  std::shared_ptr<ComPWA::FitParameter> GetMesonRadiusParameter() {
    return MesonRadius;
  }

  /// \see GetMesonRadius() const { return MesonRadius->value(); }
  void SetMesonRadius(double w) { MesonRadius->setValue(w); }

  /// Get meson radius.
  /// The meson radius is a measure of the size of the resonant state. It is
  /// used to calculate the angular momentum barrier factors.
  double GetMesonRadius() const { return MesonRadius->value(); }

  /// \see GetFormFactorType()
  void SetFormFactorType(FormFactorType t) { FFType = t; }

  /// Get form factor type.
  /// The type of formfactor that is used to calculate the angular momentum
  /// barrier factors.
  FormFactorType GetFormFactorType() { return FFType; }

  void updateParametersFrom(const ParameterList &list);
  void addUniqueParametersTo(ParameterList &list);
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<FunctionTree>
  createFunctionTree(const ParameterList &DataSample, unsigned int pos,
                     const std::string &suffix) const;

protected:
  /// Resonance spin
  ComPWA::Spin J;
  /// Orbital Angular Momentum between two daughters in Resonance decay
  ComPWA::Spin L;

  /// Masses of daughter particles
  std::pair<double, double> DaughterMasses;

  /// Names of daughter particles
  std::pair<std::string, std::string> DaughterNames;

  /// Resonance mass
  std::shared_ptr<ComPWA::FitParameter> Mass;

  /// Decay width of resonante state
  std::shared_ptr<ComPWA::FitParameter> Width;

  /// Meson radius of resonant state
  std::shared_ptr<ComPWA::FitParameter> MesonRadius;

  /// Form factor type
  FormFactorType FFType;
};

class BreitWignerStrategy : public ComPWA::Strategy {
public:
  BreitWignerStrategy(std::string namee = "")
      : ComPWA::Strategy(ParType::MCOMPLEX), name(namee) {}

  virtual const std::string to_str() const {
    return ("relativistic BreitWigner of " + name);
  }

  virtual void execute(ComPWA::ParameterList &paras,
                       std::shared_ptr<ComPWA::Parameter> &out);

protected:
  std::string name;
};

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
