//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
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
#include "Physics/HelicityFormalism/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

  class PartialDecay;
/**
 * Relativistic Breit-Wigner
 * (Breit Wigner with Blatt-Weisskopf barrier factors)
 *
 * The dynamical function implemented here is taken from PDG2014 (Eq.47-22) for
 * the one channel case.
 *
 * The three required parameters are defined in this model and can be
 * retrieved via the #getParameterList() function:
 * @param resonance_width_ width of the resonance
 * @param resonance_mass_ mass of the resonance
 * @param meson_radius_ Scale of interaction range
 *
 * The remaining required parameter is the angular momentum param J_ of the two
 * particle state.
 *
 * Additional informations from the event are extracted from the data point
 * within the #evaluate() function. For example the invariant mass and the
 * daughter masses.
 */
class RelativisticBreitWigner : public AbstractDynamicalFunction {

public:
  RelativisticBreitWigner() {};

  virtual ~RelativisticBreitWigner() {};

  std::complex<double> Evaluate(const dataPoint &point, int pos) const;

  /**! Get current normalization.  */
  virtual double GetNormalization() { return 1 / Integral(); };

  /**! Setup function tree */
  virtual std::shared_ptr<FunctionTree> SetupTree(ParameterList &sample,
                                                  ParameterList &toySample,
                                                  std::string suffix) {
    return std::shared_ptr<FunctionTree>();
  };

  // --------------------------- Set/Get functions ---------------------------

  /**
   Set decay width

   @param w Decay width
   */
  void SetWidth(std::shared_ptr<DoubleParameter> w) { _width = w; }

  /**
   Get decay width

   @return Decay width
   */
  std::shared_ptr<DoubleParameter> GetWidth() { return _width; }

  /**
   Set decay width

   @param w Decay width
   */
  void SetWidth(double w) { _width->SetValue(w); }

  /**
   Get decay width

   @return Decay width
   */
  double GetWidthValue() { return _width->GetValue(); }

  /**
   Set meson radius
   The meson radius is a measure of the size of the resonant state. It is used
   to calculate the angular momentum barrier factors.

   @param r Meson radius
   */
  void SetMesonRadius(std::shared_ptr<DoubleParameter> r) { _mesonRadius = r; }

  /**
   Get meson radius
   The meson radius is a measure of the size of the resonant state. It is used
   to calculate the angular momentum barrier factors.

   @return Meson radius
   */
  std::shared_ptr<DoubleParameter> GetMesonRadius() { return _mesonRadius; }

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
  double GetMesonRadiusValue() { return _mesonRadius->GetValue(); }

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

  /**
   Factory for RelativisticBreitWigner dynamical function

   @param pt Configuration tree
   @return Constructed object
   */
  static std::shared_ptr<RelativisticBreitWigner>
  Factory(const boost::property_tree::ptree &pt);

  std::shared_ptr<ComPWA::Physics::HelicityFormalism::PartialDecay>
  operator*(std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> wigner) {
//    auto obj =std::make_shared<ComPWA::Physics::HelicityFormalism::PartialDecay>();
//    std::shared_ptr<ComPWA::Physics::HelicityFormalism::PartialDecay>
//    partDecay(std::make_shared<PartialDecay>(this), wigner);
//    return obj;
  }

protected:
  //! Decay width of resonante state
  std::shared_ptr<DoubleParameter> _width;

  //! Meson radius of resonant state
  std::shared_ptr<DoubleParameter> _mesonRadius;

  //! Form factor type
  formFactorType _ffType;
  
private:
  //! Temporary value of meson radius (used to trigger recalculation of normalization)
  double _current_mesonRadius;
};

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_RELATIVISTICBREITWIGNER_HPP_ */
