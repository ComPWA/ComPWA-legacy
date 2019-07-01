#include "FunctionTreeEstimatorWrapper.hpp"
#include "FunctionTree.hpp"
#include "Value.hpp"

namespace ComPWA {

FunctionTreeEstimatorWrapper::FunctionTreeEstimatorWrapper(
    std::shared_ptr<FunctionTree> tree, ParameterList parameters)
    : Tree(tree), Parameters(parameters) {}

double FunctionTreeEstimatorWrapper::evaluate() {
  auto val = std::dynamic_pointer_cast<Value<double>>(Tree->parameter());
  return val->value();
}

void FunctionTreeEstimatorWrapper::updateParametersFrom(
    const std::vector<double> &params) {
  size_t pos = 0;
  for (auto p : Parameters.doubleParameters()) {
    if (!p->isFixed()) {
      p->setValue(params[pos]);
    }
    ++pos;
  }
}

std::vector<double> FunctionTreeEstimatorWrapper::getParameters() const {
  std::vector<double> params;
  for (auto p : Parameters.doubleParameters()) {
    params.push_back(p->value());
  }
  return params;
}

} /* namespace ComPWA */
