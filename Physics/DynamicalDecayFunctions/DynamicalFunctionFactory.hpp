//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef PHYSICS_DYNAMICALDECAYFUNCTIONS_DYNAMICALFUNCTIONFACTORY_HPP_
#define PHYSICS_DYNAMICALDECAYFUNCTIONS_DYNAMICALFUNCTIONFACTORY_HPP_

#include <memory>

#include "Physics/DynamicalDecayFunctions/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityAmplitude/ParticleStateDefinitions.hpp"

namespace DynamicalFunctions {

class DynamicalFunctionFactory {
  std::map<HelicityFormalism::TwoBodyDecayInformation,
      std::shared_ptr<AbstractDynamicalFunction> > dynamical_function_list_;
public:
  DynamicalFunctionFactory();
  virtual ~DynamicalFunctionFactory();

  std::shared_ptr<AbstractDynamicalFunction> generateDynamicalFunction(
      const HelicityFormalism::TwoBodyDecayInformation& state_info) const;
};

} /* namespace DynamicalFunctions */

#endif /* PHYSICS_DYNAMICALDECAYFUNCTIONS_DYNAMICALFUNCTIONFACTORY_HPP_ */
