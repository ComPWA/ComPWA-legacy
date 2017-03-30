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
#include "Physics/HelicityFormalism/AbstractDynamicalFunction.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class AmpFlatteRes : public HelicityFormalism::AbstractDynamicalFunction {
public:
  AmpFlatteRes() {};

  virtual ~AmpFlatteRes();

  virtual std::complex<double> Evaluate(const dataPoint &point, int pos) const;
  
  //! Clone function
  virtual AmpFlatteRes *Clone(std::string newName = "") const {
    auto tmp = (new AmpFlatteRes(*this));
    //if (newName != "")
    //tmp->SetName(newName);
    return tmp;
  }

  //! Check of parameters have changed and normalization has to be recalculatecd
  virtual void CheckModified() const;

  //! Print resonance parameters
  std::string to_str() const;

  /**! Get current normalization.  */
  virtual double GetNormalization() const;
  
  //! Get resonance width
  virtual double GetWidth() const {
    return std::abs(HelicityFormalism::couplingToWidth(
        _mass->GetValue(), _mass->GetValue(), _g1->GetValue(), _massA, _massB,
        (double)_spin, _mesonRadius->GetValue(), _ffType));
  }

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


  virtual std::shared_ptr<FunctionTree> GetTree(ParameterList &sample,
                                                  ParameterList &phspSample,
                                                  ParameterList &toySample,
                                                std::string suffix);

  /**
   Factory for AmpFlatteRes

   @param pt Configuration tree
   @return Constructed object
   */
  static std::shared_ptr<AmpFlatteRes>
  Factory(const boost::property_tree::ptree &pt);
  
protected:
  // Initialize masses
  void initialize();

  double _g2_massA, _g2_massB, _g3_massA, _g3_massB;
  std::shared_ptr<DoubleParameter> _g1, _g2, _g3;
  std::string _g2_idA, _g2_idB, _g3_idA, _g3_idB;

  //! Meson radius of resonant state
  std::shared_ptr<DoubleParameter> _mesonRadius;
  
  //! Form factor type
  formFactorType _ffType;
  
private:
  double _current_g3, _current_g2, _current_g1;
  
};

class FlatteStrategy : public Strategy {
public:
  FlatteStrategy(const std::string resonanceName)
      : Strategy(ParType::MCOMPLEX), name(resonanceName) {}

  virtual const std::string to_str() const {
    return ("flatte amplitude of " + name);
  }

  static std::shared_ptr<FunctionTree>
  SetupTree(std::string name, std::shared_ptr<MultiDouble> mSq,
            std::shared_ptr<DoubleParameter> mR,
            std::shared_ptr<DoubleParameter> g, double ma, double mb,
            std::shared_ptr<DoubleParameter> g2, double g2_ma, double g2_mb,
            std::shared_ptr<DoubleParameter> g3, double g3_ma, double g3_mb,
            Spin spin, std::shared_ptr<DoubleParameter> mesonRadius,
            formFactorType type);

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);

protected:
  std::string name;
};

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
