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

#ifndef PARTIALDECAY_HPP_
#define PARTIALDECAY_HPP_

#include <memory>
#include <boost/property_tree/ptree.hpp>

#include "Core/Resonance.hpp"
#include "Physics/HelicityFormalism/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"

namespace ComPWA {

std::shared_ptr<DoubleParameter>
DoubleParameterFactory(const boost::property_tree::ptree &pt) {
  auto obj = std::make_shared<DoubleParameter>();
  obj->SetValue( pt.get<double>("value") );
  return obj;
}

namespace Physics {
namespace HelicityFormalism {

class PartialDecay : ComPWA::Resonance {

public:
  PartialDecay() {};
  
  /**! Evaluate decay */
  std::complex<double> Evaluate(const dataPoint &point) const {
    std::complex<double> result =
        std::polar(_strength->GetValue(), _phase->GetValue());
    result *= _angD->Evaluate(point, 1, 2);
    result *= _dynamic->Evaluate(point);

    return result;
  };

  /**! Setup function tree */
  virtual std::shared_ptr<FunctionTree> SetupTree(ParameterList &sample,
                                                  ParameterList &toySample,
                                                  std::string suffix) {
    return std::shared_ptr<FunctionTree>();
  };

  /**
   Factory for PartialDecay

   @param pt Configuration tree
   @return Constructed object
   */
  static std::shared_ptr<PartialDecay>
  Factory(const boost::property_tree::ptree &pt) {
    auto obj = std::shared_ptr<PartialDecay>();
    obj->SetName(pt.get<string>("Resonance.<xmlattr>.name", "empty"));
    obj->SetMagnitudePar(ComPWA::DoubleParameterFactory(pt.get_child("strength")));
    obj->SetPhasePar(ComPWA::DoubleParameterFactory(pt.get_child("phase")));
    obj->SetWignerD(
        ComPWA::Physics::HelicityFormalism::AmpWignerD::Factory(pt));
    obj->SetDynamicalFunction(AbstractDynamicalFunction::Factory(pt));

    return obj;
  }

  /**
   Get WignerD function

   @return Shared_ptr<AmpWignerD>
   */
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> GetWignerD() {
    return _angD;
  }

  /**
   Set WignerD function

   @param w WignerD function
   */
  void SetWignerD(
      std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> w) {
    _angD = w;
  }

  /**
   Get dynamical function (e.g. Breit-Wigner parametrization)

   @return Shared_ptr<AbstractDynamicalFunction>
   */
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AbstractDynamicalFunction>
  SetDynamicalFunction() {
    return _dynamic;
  }

  /**
   Set dynamical function

   @param f Dynamical function
   */
  void SetDynamicalFunction(
      std::shared_ptr<
          ComPWA::Physics::HelicityFormalism::AbstractDynamicalFunction>
          f) {
    _dynamic = f;
  }

  /**
   Get strength parameter

   @return strength parameter
   */
  std::shared_ptr<ComPWA::DoubleParameter> GetMagnitudePar() { return _strength; }

  /**
   Get strength parameter

   @return strength parameter
   */
  double GetMagnitude() const { return _strength->GetValue(); }

  /**
   Set strength parameter

   @param par Strength parameter
   */
  void SetMagnitudePar(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _strength = par;
  }

  /**
   Set strength parameter

   @param par Strength parameter
   */
  void SetMagnitude(double par) { _strength->SetValue(par); }

  /**
   Get phase parameter

   @return Phase parameter
   */
  std::shared_ptr<ComPWA::DoubleParameter> GetPhasePar() { return _phase; }

  /**
   Get phase parameter

   @return Phase parameter
   */
  double GetPhase() const { return _phase->GetValue(); }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  void SetPhasePar(std::shared_ptr<ComPWA::DoubleParameter> par) { _phase = par; }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  void SetPhase(double par) { _phase->SetValue(par); }
  
  //! Get coefficient
  virtual std::complex<double> GetCoefficient() const {
//    return std::complex<double>(1,0);
  };

  //! Set prefactor
  virtual void SetPrefactor(std::complex<double> pre) { _preFactor = pre; }
  
  //! Get prefactor
  virtual std::complex<double> GetPrefactor() const { return _preFactor; }
  
  //! Implementation of interface for streaming info about the strategy
  virtual std::string to_str() const { return std::string("PartialDecay"); }

  
  
  //! Clone function
  virtual Resonance *Clone(std::string newName = "") const {
    return (new PartialDecay());
  };
  
protected:
  std::shared_ptr<ComPWA::DoubleParameter> _strength;
  std::shared_ptr<ComPWA::DoubleParameter> _phase;
  std::complex<double> _preFactor;
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> _angD;
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AbstractDynamicalFunction>
      _dynamic;

  // TODO: Operator* to construct SequentialTwoBodyDecay
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PARTIALDECAY_ */
