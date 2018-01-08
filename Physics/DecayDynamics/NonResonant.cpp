// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/DecayDynamics/NonResonant.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

std::shared_ptr<FunctionTree>
NonResonant::tree(const ParameterList &sample, int pos, std::string suffix) {

  int n = sample.mDoubleValue(0)->values().size();
  auto unitVec = MComplex("unit", n, std::complex<double>(1, 0));

  return std::make_shared<FunctionTree>("NonResonant" + suffix, unitVec);
}

} // ns::DecayDynamics
} // ns::Physics
} // ns::ComPWA
