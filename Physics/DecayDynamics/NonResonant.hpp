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
  NonResonant(std::string name = "") : AbstractDynamicalFunction(name) {
    // Set the mass parameter to make sure the pointer is set.
    SetMassParameter(
        std::shared_ptr<DoubleParameter>(new DoubleParameter("", 0.0)));
  };

  virtual ~NonResonant(){};
  
  virtual std::complex<double> evaluate(const DataPoint &p, int pos) const {
    return std::complex<double>(1.0, 0.0);
  }

  virtual std::shared_ptr<FunctionTree> tree(const ParameterList &sample,
                                                int pos, std::string suffix);

  virtual void GetParameters(ParameterList &list){};

  /// Update parameters to the values given in \p par
  virtual void updateParameters(const ParameterList &par){};
};

} // ns::DecayDynamics
} // ns::Physics
} // ns::ComPWA

#endif
