// Copyright (c) 2013, 2019 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_DYNAMICS_PRODUCTION_FORMFACTOR_HPP_
#define PHYSICS_DYNAMICS_PRODUCTION_FORMFACTOR_HPP_

#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <vector>

#include "AbstractDynamicalFunction.hpp"
#include "Core/Exceptions.hpp"
#include "Core/Functions.hpp"
#include "Core/Spin.hpp"
#include "FormFactor.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

/// The production formfactor is implemented in Blatt-Weisskopf type \f$B_{L}(q)\f$
/// Ref. Dalitz Plot Analysis, PDG 2012, where 
/// \f[
/// \frac{B^{\prime}(q_{0})}{B^{\prime}_{L}(q)} = q^{L}B^{\prime}_{L}(q, q_{0})
/// \f]
class ProductionFormFactor
    : public ComPWA::Physics::Dynamics::AbstractDynamicalFunction {

public:
  //============ CONSTRUCTION ==================
  ProductionFormFactor(std::string name,
                          std::pair<std::string, std::string> daughters,
                          std::shared_ptr<ComPWA::PartList> partL);
  virtual ~ProductionFormFactor();

  //================ EVALUATION =================

  std::complex<double> evaluate(const ComPWA::DataPoint &point,
                                unsigned int pos) const;

  /// Dynamical Blatt-Weisskopf function.
  /// \param mSq Invariant mass squared
  /// \param ma Mass of daughter particle
  /// \param mb Mass of daughter particle
  /// \param L Orbital angular momentum between two daughters a and b
  /// \param mesonRadius Meson Radius
  /// \param ffType Form factor type
  static std::complex<double>
  dynamicalFunction(double mSq, double ma, double mb, unsigned int L, 
                    double mesonRadius, FormFactorType ffType);

  //============ SET/GET =================

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
  FormFactorType GetFormFactorType() const { return FFType; }

  void updateParametersFrom(const ParameterList &list);
  void addUniqueParametersTo(ParameterList &list);
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<FunctionTree>
  createFunctionTree(const ParameterList &DataSample, unsigned int pos,
                     const std::string &suffix) const;

protected:
  /// Orbital Angular Momentum between two daughters in Resonance decay
  ComPWA::Spin L;

  /// Masses of daughter particles
  std::pair<double, double> DaughterMasses;

  /// Names of daughter particles
  std::pair<std::string, std::string> DaughterNames;

  /// Meson radius of resonant state
  std::shared_ptr<ComPWA::FitParameter> MesonRadius;

  /// Form factor type
  FormFactorType FFType;
};

class FormFactorStrategy : public ComPWA::Strategy {
public:
  FormFactorStrategy(std::string namee = "")
      : ComPWA::Strategy(ParType::MCOMPLEX), name(namee) {}

  virtual const std::string to_str() const {
    return ("production FormFactorStratey of " + name);
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
