#include "FunctionTreeEstimatorWrapper.hpp"
#include "FunctionTree.hpp"
#include "Value.hpp"

namespace ComPWA {
namespace FunctionTree {

FunctionTreeEstimatorWrapper::FunctionTreeEstimatorWrapper(
    std::shared_ptr<FunctionTree> tree, ParameterList parameters,
    ParameterList userparameters)
    : Tree(tree), Parameters(parameters), UserParameters(userparameters) {

  if (!Tree) {
    throw std::runtime_error(
        "FunctionTreeEstimatorWrapper::FunctionTreeEstimator(): "
        "FunctionTree is empty!");
  }
  Tree->parameter();
  if (!Tree->sanityCheck()) {
    throw std::runtime_error(
        "FunctionTreeEstimatorWrapper::FunctionTreeEstimatorWrapper(): Tree "
        "has structural "
        "problems. Sanity check not passed!");
  }

  // UserParameters.DeepCopy(Parameters);
  for (auto x : Parameters.doubleParameters()) {
    x->fixParameter(false);
    // IMPORTANT: we have to unfix all parameters
    // since the optimizer will take care of this
    // and we can't maintain which parameters will
    // be fixed later on. So internally we treat them
    // all as free. But of course only the parameters
    // that are not fixed will change...
    // So the function tree caching and recalculation
    // still works fine
    // When improving the FunctionTree later on
    // The Parameters should be more "dumb" and not have
    // information about a fixed or not fixed status.
  }
}

double FunctionTreeEstimatorWrapper::evaluate() {
  auto val = std::dynamic_pointer_cast<Value<double>>(Tree->parameter());
  return val->value();
}

void FunctionTreeEstimatorWrapper::updateParametersFrom(
    const std::vector<double> &params) {
  size_t pos = 0;
  for (auto p : Parameters.doubleParameters()) {
    p->setValue(params[pos]);
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

std::string FunctionTreeEstimatorWrapper::print(int level) const {
  return Tree->head()->print(level);
}

} // namespace FunctionTree
} // namespace ComPWA
