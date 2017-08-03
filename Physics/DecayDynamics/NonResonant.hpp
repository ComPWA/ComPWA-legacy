//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

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
