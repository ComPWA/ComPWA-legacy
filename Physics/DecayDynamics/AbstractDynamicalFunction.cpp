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

#include "Physics/DecayDynamics/AbstractDynamicalFunction.hpp"
#include "Tools/Integration.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

void AbstractDynamicalFunction::GetParameters(ParameterList &list) {
  /* We check of for each parameter if a parameter of the same name exists in
   *list. If so we check if both are equal and set the local parameter to the
   *parameter from the list. In this way we connect parameters that occur on
   *different positions in the amplitude.
   */
  std::shared_ptr<DoubleParameter> tmp;
  auto mass = GetMassParameter();
  try { // catch BadParameter
    tmp = list.GetDoubleParameter(mass->GetName());
    try { //catch and throw std::runtime_error due to failed parameter comparisson
      if (*tmp == *mass)
        SetMassParameter(tmp);
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.AddParameter(mass);
  }
}

} /* namespace DecayDynamics */
} /* namespace Physics */
} /* namespace ComPWA */
