// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/DecayDynamics/AbstractDynamicalFunction.hpp"
#include "Tools/Integration.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

void AbstractDynamicalFunction::GetParameters(ParameterList &list) {
  // We check of for each parameter if a parameter of the same name exists in
  // list. If so we check if both are equal and set the local parameter to the
  // parameter from the list. In this way we connect parameters that occur on
  // different positions in the amplitude.
  std::shared_ptr<FitParameter> tmp;
  auto mass = GetMassParameter();
  try { // catch BadParameter
    tmp = FindParameter(mass->name(), list);
    try { //catch and throw std::runtime_error due to failed parameter comparisson
      if (*tmp == *mass)
        SetMassParameter(tmp);
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.addParameter(mass);
  }
}

} // namespace DecayDynamics
} // namespace Physics
} // namespace ComPWA
