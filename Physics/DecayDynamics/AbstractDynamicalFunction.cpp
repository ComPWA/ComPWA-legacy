// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/DecayDynamics/AbstractDynamicalFunction.hpp"
#include "Tools/Integration.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

void AbstractDynamicalFunction::parameters(ParameterList &list) {
  // We check of for each parameter if a parameter of the same name exists in
  // list. If so we check if both are equal and set the local parameter to the
  // parameter from the list. In this way we connect parameters that occur on
  // different positions in the amplitude.
 Mass = list.addUniqueParameter(Mass);
}

} // namespace DecayDynamics
} // namespace Physics
} // namespace ComPWA
