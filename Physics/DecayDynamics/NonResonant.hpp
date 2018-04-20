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
  //============ CONSTRUCTION ==================
  NonResonant(std::string name = "") : AbstractDynamicalFunction(name) {
    // Set the mass parameter to make sure the pointer is set.
    SetMassParameter(
        std::shared_ptr<FitParameter>(new FitParameter("", 0.0)));
  };

  virtual ~NonResonant(){};

  //======= INTEGRATION/NORMALIZATION ===========
  /// Check of parameters have changed and normalization has to be
  /// recalculatecd. Since this function does not have any parameters false is
  /// returned in all cases.
  virtual bool isModified() const { return false; }

  /// Label as modified/unmodified. Empty for this function.
  virtual void setModified(bool b) {};

  //================ EVALUATION =================
  
  virtual std::complex<double> evaluate(const DataPoint &p, int pos) const {
    return std::complex<double>(1.0, 0.0);
  }

  virtual void parameters(ParameterList &list){};

  /// Update parameters to the values given in \p par
  virtual void updateParameters(const ParameterList &par){};

  virtual std::shared_ptr<FunctionTree> tree(const ParameterList &sample,
                                             int pos, std::string suffix);
};

} // ns::DecayDynamics
} // ns::Physics
} // ns::ComPWA

#endif
