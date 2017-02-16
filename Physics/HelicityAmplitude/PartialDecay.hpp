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

#ifndef PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDE_HPP_
#define PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDE_HPP_

#inlcude "Core/Resonance.hpp"
#include "Physics/DynamicalDecayFunctions/AbstractDynamicalFunction.hpp"
#include "Physics/AmpitudeSum/AmpWignerD.hpp"

namespace ComPWA {
std::shared_ptr<DoubleParameter>
DoubleParameterFactory(boost::property_tree::ptree &pt) {
  auto obj = std::make_shared<DoubleParameter>();
    obj->SetValue( pt.second.get_child("value");
}
namespace Physics {

namespace HelicityFormalism {

class PartialDecay : Resonance {

public:
  /**! Evaluate decay */
  std::complex<double> Evaluate(dataPoint &point) {
    std::complex<double> result =
        std::polar(strength_->GetValue(), phase_->GetValue());
    result *= angular_part_->Evaluate(point);
    result *= dynamical_part_->Evaluate(point);

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
  Factory(boost::property_tree::ptree &pt) {
    auto obj = std::make_shared<PartialDecay>();
    obj->SetName(pt.get<string>("Resonance.<xmlattr>.name", "empty"));
    obj->SetStrength(DoubleParameterFactory(pt.get_child("strength")));
    obj->SetPhase(DoubleParameterFactory(pt.get_child("phase")));
    obj->SetWignerD(AmpWignerD::Factory(pt));
    obj->SetDynamicalFunction(AbstractDynamicalFunction::Factory(pt));

    return obj;
  }

  /**
   Get WignerD function

   @return Shared_ptr<AmpWignerD>
   */
  std::shared_ptr<ComPWA::Physics::AmpSumIntensity::AmpWignerD> GetWignerD() {
    return _angD;
  }

  /**
   Set WignerD function

   @param w WignerD function
   */
  void
  SetWignerD(std::shared_ptr<ComPWA::Physics::AmpSumIntensity::AmpWignerD> w) {
    _angD = w;
  }

  /**
   Get dynamical function (e.g. Breit-Wigner parametrization)

   @return Shared_ptr<AbstractDynamicalFunction>
   */
  std::shared_ptr<
      ComPWA::Physics::DynamicalFunctions::AbstractDynamicalFunction>
  SetDynamicalFunction() {
    return _dynamic;
  }

  /**
   Set dynamical function

   @param f Dynamical function
   */
  void SetDynamicalFunction(
      std::shared_ptr<
          ComPWA::Physics::DynamicalFunctions::AbstractDynamicalFunction>
          f) {
    _dynamic = f;
  }

  /**
   Get strength parameter

   @return strength parameter
   */
  std::shared_ptr<ComPWA::Physics::AmpSumIntensity::AmpWignerD> GetStrength() {
    return _strength;
  }

  /**
   Get strength parameter

   @return strength parameter
   */
  double GetStrength() { return _strength->GetValue(); }

  /**
   Set strength parameter

   @param par Strength parameter
   */
  void SetStrength(std::shared_ptr<ComPWA::DoubleParameter> par) {
    _strength = par;
  }

  /**
   Set strength parameter

   @param par Strength parameter
   */
  void SetStrength(double par) { _strength->SetValue(par); }

  /**
   Get phase parameter

   @return Phase parameter
   */
  std::shared_ptr<ComPWA::Physics::AmpSumIntensity::AmpWignerD> GetPhase() {
    return _phase;
  }

  /**
   Get phase parameter

   @return Phase parameter
   */
  double GetPhase() { return _phase->GetValue(); }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  void SetPhase(std::shared_ptr<ComPWA::DoubleParameter> par) { _phase = par; }

  /**
   Set phase parameter

   @param par Phase parameter
   */
  void SetPhase(double par) { _phase->SetValue(par); }

protected:
  std::shared_ptr<ComPWA::DoubleParameter> _strength;
  std::shared_ptr<ComPWA::DoubleParameter> _phase;
  std::shared_ptr<ComPWA::Physics::AmpSumIntensity::AmpWignerD> _angD;
  std::shared_ptr<
      ComPWA::Physics::DynamicalFunctions::AbstractDynamicalFunction>
      _dynamic;

  // TODO: Operator* to construct SequentialTwoBodyDecay
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_TOPOLOGYAMPLITUDE_HPP_ */

/
