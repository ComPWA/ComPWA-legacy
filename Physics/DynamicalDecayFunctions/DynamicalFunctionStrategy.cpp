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
    std::shared_ptr<AbstractDynamicalFunction> abs_dyn_func,
    unsigned int storage_index) :
    Strategy(ParType::MCOMPLEX), dynamical_function_(abs_dyn_func), storage_index_(
        storage_index) {
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
    values.reserve(DataPointStorage::Instance().getNumberOfEvents(storage_index_));
    for (unsigned int i = 0;
        i < DataPointStorage::Instance().getNumberOfEvents(storage_index_); ++i) {
      values.push_back(
          dynamical_function_->evaluate(storage_index_, i,
              paras.GetMultiUnsignedInteger(0)->GetValue(0)));
    }
    out = std::shared_ptr<MultiComplex>(
        new MultiComplex(out->GetName(), values));
  }
  return true;
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
