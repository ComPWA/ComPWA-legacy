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

#ifndef PHYSICS_DYNAMICALDECAYFUNCTIONS_TWOBODYDECAY_TOPNODECONSTANTVALUE_HPP_
#define PHYSICS_DYNAMICALDECAYFUNCTIONS_TWOBODYDECAY_TOPNODECONSTANTVALUE_HPP_

#include "Physics/HelicityFormalism/AbstractDynamicalFunction.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class NonResonant : public AbstractDynamicalFunction {

public:
  NonResonant() {};
  virtual ~NonResonant() {};

  virtual std::complex<double> Evaluate(const dataPoint &p) const {
    return EvaluateNoNorm(0.0) * GetNormalization();
  }
  virtual std::complex<double> EvaluateNoNorm(double mSq) const {
    return std::complex<double>(1.0, 0.0);
  }

  /**! Get current normalization.  */
  virtual double GetNormalization() const { return (1 / Integral()); };

  virtual std::shared_ptr<FunctionTree> GetTree(ParameterList &sample,
                                        ParameterList &toySample,
                                        std::string suffix);

  virtual void GetParameters(ParameterList &list) {};
  
protected:
  virtual double Integral() const {
    return Kinematics::Instance()->GetPhspVolume();
  }
};

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_DYNAMICALDECAYFUNCTIONS_TWOBODYDECAY_TOPNODECONSTANTVALUE_HPP_ \
          */
