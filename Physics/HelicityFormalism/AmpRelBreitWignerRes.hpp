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
//handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef AMP_REL_BREIT_WIGNER_RES
#define AMP_REL_BREIT_WIGNER_RES

#include <vector>

#include "Core/Functions.hpp"
#include "Core/Exceptions.hpp"
#include "Core/FunctionTree.hpp"
#include "Physics/HelicityFormalism/AbstractDynamicalFunction.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class AmpRelBreitWignerRes : public AbstractDynamicalFunction {
public:
  AmpRelBreitWignerRes() {};

  virtual ~AmpRelBreitWignerRes();

  virtual std::complex<double> Evaluate(const dataPoint &point, int pos) const;
  
  //! Clone function
  virtual AmpRelBreitWignerRes *Clone(std::string newName = "") const {
    auto tmp = (new AmpRelBreitWignerRes(*this));
//    if (newName != "")
//      tmp->SetName(newName);
    return tmp;
  }

  //! Trigger recalculation of normalization
  virtual void CheckModified() const;

  std::string to_str() const;

  //! Calculation integral |dynamical amplitude|^2
  virtual double GetNormalization() const;

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
  
   // --------------------------------------------------------------------
  
  /** Breit-Wigner function
   *
   * The dynamical function implemented here is taken from PDG2014 (Eq.47-22)
   * for
   * the one channel case.
   * @J is angular momentum between A&B. In case A&B have spin 0 this is the
   * spin for the resonace.
   *
   * @param mSq Invariant mass
   * @param mR Resonance mass
   * @param ma Mass particle A
   * @param mb Mass particle B
   * @param width  Width of resonance
   * @param J Angular momentum between A&B
   * @param mesonRadius Scale of interaction range
   * @return
   */
  static std::complex<double>
  dynamicalFunction(double mSq, double mR, double ma, double mb, double width,
                    unsigned int J, double mesonRadius,
                    formFactorType ffType = formFactorType::BlattWeisskopf);


  virtual std::shared_ptr<FunctionTree> GetTree(ParameterList &sample,
                                                  ParameterList &phspSample,
                                                  ParameterList &toySample,
                                                  std::string suffix);

protected:
  //! Resonance width
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

//class BreitWignerStrategy : public Strategy {
//public:
//  BreitWignerStrategy(const std::string resonanceName)
//      : Strategy(ParType::MCOMPLEX), name(resonanceName) {}
//
//  virtual const std::string to_str() const {
//    return ("relativistic BreitWigner of " + name);
//  }
//
//  virtual bool execute(ParameterList &paras,
//                       std::shared_ptr<AbsParameter> &out);
//
//protected:
//  std::string name;
//};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
