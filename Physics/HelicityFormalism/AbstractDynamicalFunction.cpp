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

#include "Physics/HelicityFormalism/AbstractDynamicalFunction.hpp"
#include "Tools/Integration.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

//! Integral
double AbstractDynamicalFunction::Integral() const {
  if (!_phspSample->size()) {
    LOG(debug)
        << "CoherentIntensity::Integral() | Integral can not be calculated "
           "since no phsp sample is set. Set a sample using "
           "SetPhspSamples(phspSample, toySample)!";
    return 1.0;
  }

  double sumIntens = 0;
  for (auto i : *_phspSample.get())
    sumIntens += std::norm(EvaluateNoNorm(i.GetValue(_dataPos)));
  
  double phspVol = Kinematics::Instance()->GetPhspVolume();
  double integral = (sumIntens * phspVol / _phspSample->size());
  LOG(trace) << "AbstractDynamicalFunction::Integral() | Integral is " << integral
             << ".";

  return integral;

//Numeric integration in (1-dim) invariant mass range
//  double phspVol = Kinematics::Instance()->GetPhspVolume();
//  double integral = (sumIntens * phspVol / _phspSample->size());
//  LOG(trace) << "CoherentIntensity::Integral() | Integral is " << integral
//             << " and the maximum value of intensity is " << maxVal << ".";
//
//  const_cast<double &>(_maxIntens) = maxVal;
//  return integral;
//  
//  auto intAlg = ComPWA::Tools::IntegralByQuadrature<AbstractDynamicalFunction>(
//      *this, _limits);
//  double integral = intAlg.Integral();
//  LOG(trace) << "AbstractDynamicalFunction::Integral() | Integral is "
//             << integral << ".";
//  return std::sqrt(integral);
}

void AbstractDynamicalFunction::GetParameters(ParameterList &list) {
  /* We check of for each parameter if a parameter of the same name exists in
   *list. If so we check if both are equal and set the local parameter to the
   *parameter from the list. In this way we connect parameters that occur on
   *different positions in the amplitude.
   */
  std::shared_ptr<DoubleParameter> tmp;
  auto mass = GetMass();
  try { // catch BadParameter
    tmp = list.GetDoubleParameter(mass->GetName());
    try { //catch and throw std::runtime_error due to failed parameter comparisson
      if (*tmp == *mass)
        SetMass(tmp);
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.AddParameter(mass);
  }
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
