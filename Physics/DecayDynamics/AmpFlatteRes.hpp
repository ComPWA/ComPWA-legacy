//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding correct g1s
//-------------------------------------------------------------------------------
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
  AmpFlatteRes() : AbstractDynamicalFunction(){};

  virtual ~AmpFlatteRes();

  //! Clone function
  virtual AmpFlatteRes *Clone(std::string newName = "") const {
    auto tmp = (new AmpFlatteRes(*this));
    // if (newName != "")
    // tmp->SetName(newName);
    return tmp;
  }

  /**
   Factory for AmpFlatteRes

   @param pt Configuration tree
   @return Constructed object
   */
  static std::shared_ptr<AbstractDynamicalFunction>
  Factory(const boost::property_tree::ptree &pt);

  //======= INTEGRATION/NORMALIZATION ===========

  //! Check of parameters have changed and normalization has to be recalculatecd
  virtual bool CheckModified() const;

  //================ EVALUATION =================

  virtual std::complex<double> Evaluate(const dataPoint &point, int pos) const;

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

  virtual void GetParameters(ParameterList &list);

  //! Fill vector with parameters
  virtual void GetParametersFast(std::vector<double> &list) const {
    AbstractDynamicalFunction::GetParametersFast(list);
    for (auto i : _g)
      list.push_back(i.GetValue());
    list.push_back(GetMesonRadius());
  }

  /**
  Set meson radius
  The meson radius is a measure of the size of the resonant state. It is used
  to calculate the angular momentum barrier factors.

  @param r Meson radius
  */
  void SetMesonRadiusParameter(std::shared_ptr<DoubleParameter> r) {
    _mesonRadius = r;
  }

  /**
   Get meson radius
   The meson radius is a measure of the size of the resonant state. It is used
   to calculate the angular momentum barrier factors.

   @return Meson radius
   */
  std::shared_ptr<DoubleParameter> GetMesonRadiusParameter() {
    return _mesonRadius;
  }

  /**
   Set meson radius
   The meson radius is a measure of the size of the resonant state. It is used
   to calculate the angular momentum barrier factors.

   @param r Meson radius
   */
  void SetMesonRadius(double w) { _mesonRadius->SetValue(w); }

  /**
   Get meson radius
   The meson radius is a measure of the size of the resonant state. It is used
   to calculate the angular momentum barrier factors.

   @return Meson radius
   */
  double GetMesonRadius() const { return _mesonRadius->GetValue(); }

  /**
   Set form factor type
   The type of formfactor that is used to calculate the angular momentum barrier
   factors

   @param t From factor type
   */
  void SetFormFactorType(formFactorType t) { _ffType = t; }

  /**
   Get form factor type
   The type of formfactor that is used to calculate the angular momentum barrier
   factors

   @return From factor type
   */
  formFactorType GetFormFactorType() { return _ffType; }

  /*! Set coupling parameter to signal channel.
  */
  void SetCoupling(Coupling g1, Coupling g2 = Coupling(),
                   Coupling g3 = Coupling()) {
    _g = std::vector<Coupling>{g1, g2, g3};
  }

  /*! Get coupling parameter to channel @channel.
   */
  Coupling GetCoupling(int channel) { return _g.at(channel); }

  /*! Get coupling parameter.
   */
  std::vector<Coupling> GetCouplings(int i) const { return _g; }

  /*! Set coupling parameters.
   */
  void SetCouplings(std::vector<Coupling> vC) {
    if (vC.size() != 2 && vC.size() != 3)
      throw std::runtime_error(
          "AmpFlatteRes::SetCouplings() | Vector with "
          "couplings has a wrong size. We expect either 2 or 3 couplings.");

    _g = vC;

    if (_g.size() == 2)
      _g.push_back(Coupling());
    // Check if one of the  coupling match the final state (_daughterMasses)
    auto mm = GetDecayMasses();
    if (mm == std::pair<double, double>(-999, -999))
      LOG(info)
          << "AmpFlatteRes::SetCouplings() | Masses of decay products not set. "
             " Can not determine if correct couplings were set.";

    bool ok = false;
    for (auto i : _g) {
      if (i.GetMassA() == mm.first && i.GetMassB() == mm.second)
        ok = true;
      if (i.GetMassB() == mm.first && i.GetMassA() == mm.second)
        ok = true;
    }
    if (!ok)
      throw std::runtime_error("AmpFlatteRes::SetCouplings() | No couplings "
                               "for the current decay particles set!");
  }

  //=========== FUNCTIONTREE =================

  virtual std::shared_ptr<FunctionTree> GetTree(const ParameterList &sample,
                                                int pos, std::string suffix);

protected:
  // Initialize masses
  void initialize();

  //! Meson radius of resonant state
  std::shared_ptr<DoubleParameter> _mesonRadius;

  std::vector<Coupling> _g;

  //! Form factor type
  formFactorType _ffType;

private:
  double _current_g, _current_gHidden, _current_gHidden2;
};

class FlatteStrategy : public Strategy {
public:
  FlatteStrategy(const std::string resonanceName)
      : Strategy(ParType::MCOMPLEX), name(resonanceName) {}

  virtual const std::string to_str() const {
    return ("flatte amplitude of " + name);
  }

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);

protected:
  std::string name;
};

} /* namespace DecayDynamics */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
