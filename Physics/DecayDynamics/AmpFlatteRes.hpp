// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef AMP_FLATTE3CH_RES
#define AMP_FLATTE3CH_RES

#include <vector>
#include <cmath>
#include "Physics/DecayDynamics/AbstractDynamicalFunction.hpp"
#include "Physics/DecayDynamics/FormFactor.hpp"
#include "Physics/DecayDynamics/Coupling.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

class AmpFlatteRes : public DecayDynamics::AbstractDynamicalFunction {

public:
  //============ CONSTRUCTION ==================
  AmpFlatteRes()
      : AbstractDynamicalFunction(), FormFactorType(noFormFactor), Current_g(0.0),
        Current_gHidden(0.0), Current_gHidden2(0.0){};

  AmpFlatteRes(std::string name, std::pair<std::string,std::string> daughters,
               std::shared_ptr<ComPWA::PartList> partL);

  virtual ~AmpFlatteRes();

  virtual AmpFlatteRes *Clone(std::string newName = "") const {
    auto tmp = (new AmpFlatteRes(*this));
    return tmp;
  }

  //======= INTEGRATION/NORMALIZATION ===========

  //! Check of parameters have changed and normalization has to be recalculatecd
  virtual bool isModified() const;

  //================ EVALUATION =================

  virtual std::complex<double> evaluate(const DataPoint &point, int pos) const;

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
   * @param J resonance spin
   * @param mesonRadius 1/interaction length (needed for barrier factors)
   * @param ffType formfactor type
   * @return
   */
  static std::complex<double>
  dynamicalFunction(double mSq, double mR, double massA1, double massA2,
                    double gA, double massB1, double massB2, double gB,
                    double massC1, double massC2, double gC, unsigned int J,
                    double mesonRadius, formFactorType ffType);

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

  void SetMesonRadiusParameter(std::shared_ptr<FitParameter> r) {
    MesonRadius = r;
  }

  std::shared_ptr<FitParameter> GetMesonRadiusParameter() {
    return MesonRadius;
  }

  void SetMesonRadius(double w) { MesonRadius->setValue(w); }

  double GetMesonRadius() const { return MesonRadius->value(); }

  void SetFormFactorType(formFactorType t) { FormFactorType = t; }

  formFactorType GetFormFactorType() { return FormFactorType; }

  /// Set coupling parameter to signal channel and up to two more hidden
  /// channels.
  void SetCoupling(Coupling g1, Coupling g2 = Coupling(0.0, 0.0, 0.0),
                   Coupling g3 = Coupling(0.0, 0.0, 0.0)) {
    Couplings = std::vector<Coupling>{g1, g2, g3};
  }

  Coupling GetCoupling(int channel) { return Couplings.at(channel); }

  std::vector<Coupling> GetCouplings(int i) const { return Couplings; }

  void SetCouplings(std::vector<Coupling> vC);

  virtual void parameters(ParameterList &list);

  /// Fill vector with parameters
  virtual void parametersFast(std::vector<double> &list) const {
    AbstractDynamicalFunction::parametersFast(list);
    for (auto i : Couplings)
      list.push_back(i.value());
    list.push_back(GetMesonRadius());
  }

  /// Update parameters to the values given in \p par
  virtual void updateParameters(const ParameterList &par);

  //=========== FUNCTIONTREE =================

  virtual std::shared_ptr<FunctionTree> tree(const ParameterList &sample,
                                                int pos, std::string suffix);

protected:
  /// Meson radius of resonant state
  std::shared_ptr<FitParameter> MesonRadius;

  /// Coupling parameters and final state masses for multiple channels
  std::vector<Coupling> Couplings;

  /// Form factor type
  formFactorType FormFactorType;

private:
  double Current_g, Current_gHidden, Current_gHidden2;
};

class FlatteStrategy : public Strategy {
public:
  FlatteStrategy(const std::string resonanceName)
      : Strategy(ParType::MCOMPLEX), name(resonanceName) {}

  virtual const std::string to_str() const {
    return ("flatte amplitude of " + name);
  }

  virtual void execute(ParameterList &paras,
                       std::shared_ptr<Parameter> &out);

protected:
  std::string name;
};

} // ns::DecayDynamics
} // ns::Physics
} // ns::ComPWA

#endif
