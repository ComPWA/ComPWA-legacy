// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_DYNAMICS_FLATTE_HPP_
#define COMPWA_PHYSICS_DYNAMICS_FLATTE_HPP_

#include <cmath>
#include <vector>

#include "AbstractDynamicalFunction.hpp"
#include "Core/Properties.hpp"
#include "Coupling.hpp"
#include "FormFactor.hpp"

namespace ComPWA {
namespace Physics {
namespace Dynamics {

class Flatte : public AbstractDynamicalFunction {

public:
  //============ CONSTRUCTION ==================
  Flatte() : AbstractDynamicalFunction(), FFType(noFormFactor){};

  Flatte(std::string name, std::pair<std::string, std::string> daughters,
         std::shared_ptr<ComPWA::PartList> partL);

  virtual ~Flatte();

  boost::property_tree::ptree save() const {
    return boost::property_tree::ptree();
  }

  //================ EVALUATION =================

  std::complex<double> evaluate(const ComPWA::DataPoint &point,
                                unsigned int pos) const;

  /** Dynamical function for two coupled channel approach
   *
   * @param mSq center-of-mass energy^2 (=s)
   * @param mR mass of resonances
   * @param massA1 mass of first particle of signal channel
   * @param massA2 mass of second particle of signal channel
   * @param gA coupling constant for signal channel
   * @param massB1 mass of first particle of second channel
   * @param massB2 mass of second particle of second channel
   * @param gB coupling constant for second channel
   * @param massC1 mass of first particle of third channel
   * @param massC2 mass of third particle of third channel
   * @param gC coupling constant for third channel
   * @param L Orbital angular momentum between two daughters a and b
   * @param mesonRadius 1/interaction length (needed for barrier factors)
   * @param ffType formfactor type
   * @return
   */
  static std::complex<double>
  dynamicalFunction(double mSq, double mR, double massA1, double massA2,
                    double gA, double massB1, double massB2, double gB,
                    double massC1, double massC2, double gC, unsigned int J,
                    double mesonRadius, FormFactorType ffType);

  /** Dynamical function for two coupled channel approach
   *
   * @param mSq center-of-mass energy^2 (=s)
   * @param mR mass of resonances
   * @param gA coupling constant for signal channel
   * @param termA Coupling term to signal channel
   * @param termB Coupling term to second channel
   * @param termC Coupling term to third channel (optional)
   * @return
   */
  static std::complex<double>
  dynamicalFunction(double mSq, double mR, double gA,
                    std::complex<double> termA, std::complex<double> termB,
                    std::complex<double> termC = std::complex<double>(0, 0));

  //============ SET/GET =================

  void SetOrbitalAngularMomentum(const ComPWA::Spin &L_) { L = L_; }

  void SetMesonRadiusParameter(std::shared_ptr<FitParameter> r) {
    MesonRadius = r;
  }

  std::shared_ptr<FitParameter> GetMesonRadiusParameter() {
    return MesonRadius;
  }

  void SetMesonRadius(double w) { MesonRadius->setValue(w); }

  double GetMesonRadius() const { return MesonRadius->value(); }

  void SetFormFactorType(FormFactorType t) { FFType = t; }

  FormFactorType GetFormFactorType() { return FFType; }

  /// Set coupling parameter to signal channel and up to two more hidden
  /// channels.
  void SetCoupling(Coupling g1, Coupling g2 = Coupling(0.0, 0.0, 0.0),
                   Coupling g3 = Coupling(0.0, 0.0, 0.0)) {
    Couplings = std::vector<Coupling>{g1, g2, g3};
  }

  Coupling GetCoupling(int channel) { return Couplings.at(channel); }

  std::vector<Coupling> GetCouplings(int i) const { return Couplings; }

  void SetCouplings(std::vector<Coupling> vC);

  void updateParametersFrom(const ParameterList &list);
  void addUniqueParametersTo(ParameterList &list);
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  //=========== FUNCTIONTREE =================

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

  /// Resonance mass
  std::shared_ptr<ComPWA::FitParameter> Mass;

  /// Meson radius of resonant state
  std::shared_ptr<FitParameter> MesonRadius;

  /// Coupling parameters and final state masses for multiple channels
  std::vector<Coupling> Couplings;

  /// Form factor type
  FormFactorType FFType;
};

class FlatteStrategy : public Strategy {
public:
  FlatteStrategy(const std::string resonanceName)
      : Strategy(ParType::MCOMPLEX), name(resonanceName) {}

  virtual const std::string to_str() const {
    return ("flatte amplitude of " + name);
  }

  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out);

protected:
  std::string name;
};

} // namespace Dynamics
} // namespace Physics
} // namespace ComPWA

#endif
