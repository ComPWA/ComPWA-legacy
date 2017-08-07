//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//    Mathias Michel - initial API and implementation
//		Peter Weidenkaff - correct nominator, using dataPoint for data
// handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef PHYSICS_HELICITYAMPLITUDE_RELATIVISTICBREITWIGNER_HPP_
#define PHYSICS_HELICITYAMPLITUDE_RELATIVISTICBREITWIGNER_HPP_

#include <vector>
#include <memory>
#include <boost/property_tree/ptree.hpp>

#include "Core/Spin.hpp"
#include "Core/Functions.hpp"
#include "Core/Exceptions.hpp"
#include "Physics/DecayDynamics/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

class PartialDecay;
/**
 * Relativistic Breit-Wigner
 * (Breit Wigner with Blatt-Weisskopf barrier factors)
 *
 * The dynamical function implemented here is taken from PDG2014 (Eq.47-22) for
 * the one channel case.
 *
 * The three required parameters are defined in this model and can be
 * retrieved via the #GetParameters() function:
 * @param resonance_width_ Width of the resonance
 * @param resonance_mass_ Mass of the resonance
 * @param meson_radius_ Scale of interaction range
 */
class RelativisticBreitWigner : public AbstractDynamicalFunction {

public:
  //============ CONSTRUCTION ==================

  RelativisticBreitWigner(){};

  virtual ~RelativisticBreitWigner(){};

  /**
   Factory for RelativisticBreitWigner

   @param pt Configuration tree
   @return Constructed object
   */
  static std::shared_ptr<AbstractDynamicalFunction>
  Factory(const boost::property_tree::ptree &pt);

  //======= INTEGRATION/NORMALIZATION ===========

  //! Check of parameters have changed and normalization has to be recalculatecd
  virtual bool CheckModified() const;

  //================ EVALUATION =================
  std::complex<double> Evaluate(const dataPoint &point, int pos) const;

  /**
   Dynamical Breit-Wigner function

   @param mSq Invariant mass squared
   @param mR Mass of the resonant state
   @param ma Mass of daughter particle
   @param mb Mass of daughter particle
   @param width Decay width
   @param J Spin
   @param mesonRadius Meson Radius
   @param ffType Form factor type
   @return Amplitude value
   */
  static std::complex<double>
  dynamicalFunction(double mSq, double mR, double ma, double mb, double width,
                    unsigned int J, double mesonRadius, formFactorType ffType);

  //============ SET/GET =================

  /**
   Set decay width

   @param w Decay width
   */
  void SetWidthParameter(std::shared_ptr<DoubleParameter> w) { _width = w; }

  /**
   Get decay width

   @return Decay width
   */
  std::shared_ptr<DoubleParameter> GetWidthParameter() { return _width; }

  /**
   Set decay width

   @param w Decay width
   */
  void SetWidth(double w) { _width->SetValue(w); }

  /**
   Get decay width

   @return Decay width
   */
  double GetWidth() const { return _width->GetValue(); }

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

  virtual void GetParameters(ParameterList &list);

  //! Fill vector with parameters
  virtual void GetParametersFast(std::vector<double> &list) const {
    AbstractDynamicalFunction::GetParametersFast(list);
    list.push_back(GetWidth());
    list.push_back(GetMesonRadius());
  }

  //=========== FUNCTIONTREE =================

  //! Check of tree is available
  virtual bool HasTree() const { return true; }

  /**! Setup function tree */
  virtual std::shared_ptr<FunctionTree>
  GetTree(const ParameterList &sample, int pos, std::string suffix = "");

protected:
  //! Decay width of resonante state
  std::shared_ptr<DoubleParameter> _width;

  //! Meson radius of resonant state
  std::shared_ptr<DoubleParameter> _mesonRadius;

  //! Form factor type
  formFactorType _ffType;

private:
  //! Temporary values (used to trigger recalculation of normalization)
  double _current_mesonRadius;
  double _current_width;
};

class BreitWignerStrategy : public Strategy {
public:
  BreitWignerStrategy(const std::string resonanceName)
      : Strategy(ParType::MCOMPLEX), name(resonanceName) {}

  virtual const std::string to_str() const {
    return ("relativistic BreitWigner of " + name);
  }

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);

protected:
  std::string name;
};

} /* namespace DecayDynamics */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_RELATIVISTICBREITWIGNER_HPP_ */
