/*
 * DynamicalFunctionStrategy.cpp
 *
 *  Created on: Aug 9, 2016
 *      Author: steve
 */

#include "Core/DataPointStorage.hpp"
#include "Physics/DynamicalDecayFunctions/DynamicalFunctionStrategy.hpp"

namespace ComPWA {
namespace Physics {
namespace DynamicalFunctions {

DynamicalFunctionStrategy::DynamicalFunctionStrategy(
    std::shared_ptr<AbstractDynamicalFunction> abs_dyn_func) :
    Strategy(ParType::MCOMPLEX), dynamical_function_(abs_dyn_func) {
}

DynamicalFunctionStrategy::~DynamicalFunctionStrategy() {
}

//! Pure Virtual interface for streaming info about the strategy
const std::string DynamicalFunctionStrategy::to_str() const {
  return "DynFunc";
}

//! Pure Virtual interface for executing a strategy
bool DynamicalFunctionStrategy::execute(ParameterList& paras,
    std::shared_ptr<AbsParameter>& out) {

  if (ParType::MCOMPLEX) {
    std::vector<std::complex<double> > values;
    values.reserve(DataPointStorage::Instance().getNumberOfEvents());
    for (unsigned int i = 0;
        i < DataPointStorage::Instance().getNumberOfEvents(); ++i) {
      std::complex<double> value(0.0);
      auto const& evaluation_list =
          paras.GetMultiUnsignedIntegers()[0]->GetValues();
      //should be just one entry in that list
      for (unsigned int j = 0; j < evaluation_list.size(); ++j) {
        value += dynamical_function_->evaluate(i, evaluation_list[j]);
      }
      values.push_back(value);
    }
    out = std::shared_ptr<MultiComplex>(new MultiComplex(out->GetName(), values));
  }
  return true;
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
