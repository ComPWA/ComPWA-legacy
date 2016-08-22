/*
 * DynamicalFunctionStrategy.hpp
 *
 *  Created on: Aug 9, 2016
 *      Author: steve
 */

#ifndef PHYSICS_DYNAMICALDECAYFUNCTIONS_DYNAMICALFUNCTIONSTRATEGY_HPP_
#define PHYSICS_DYNAMICALDECAYFUNCTIONS_DYNAMICALFUNCTIONSTRATEGY_HPP_

#include "Core/Functions.hpp"
#include "Physics/DynamicalDecayFunctions/AbstractDynamicalFunction.hpp"

namespace ComPWA {
namespace Physics {
namespace DynamicalFunctions {

class DynamicalFunctionStrategy: public Strategy {
  std::shared_ptr<AbstractDynamicalFunction> dynamical_function_;

public:
  DynamicalFunctionStrategy(
      std::shared_ptr<AbstractDynamicalFunction> abs_dyn_func);
  virtual ~DynamicalFunctionStrategy();

  //! Pure Virtual interface for streaming info about the strategy
  virtual const std::string to_str() const;

  //! Pure Virtual interface for executing a strategy
  virtual bool execute(ParameterList& paras,
      std::shared_ptr<AbsParameter>& out);

};

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_DYNAMICALDECAYFUNCTIONS_DYNAMICALFUNCTIONSTRATEGY_HPP_ */
