// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_HELICITYFORMALISM_NONRESONANT
#define PHYSICS_HELICITYFORMALISM_NONRESONANT

#include "Physics/DecayDynamics/AbstractDynamicalFunction.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

class NonResonant : public AbstractDynamicalFunction {

public:
  NonResonant() {
    // Set the mass parameter to make shure the pointer is set.
    SetMassParameter(
        std::shared_ptr<DoubleParameter>(new DoubleParameter("", 0.0)));
  };

  virtual ~NonResonant(){};

  virtual std::complex<double> Evaluate(const dataPoint &p, int pos) const {
    return std::complex<double>(1.0, 0.0);
  }

  virtual std::shared_ptr<FunctionTree> GetTree(const ParameterList &sample,
                                                int pos, std::string suffix);

  virtual void GetParameters(ParameterList &list){};
};

} /* namespace DecayDynamics */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYFORMALISM_NONRESONANT */
