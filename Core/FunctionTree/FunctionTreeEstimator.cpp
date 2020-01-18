#include "FunctionTreeEstimator.hpp"

#include "TreeNode.hpp"
#include "Value.hpp"

namespace ComPWA {
namespace FunctionTree {

FunctionTreeEstimator::FunctionTreeEstimator(std::shared_ptr<TreeNode> tree,
                                             ParameterList parameters)
    : Tree(tree), Parameters(parameters) {

  if (!Tree) {
    throw std::runtime_error("FunctionTreeEstimator::FunctionTreeEstimator(): "
                             "FunctionTree is empty!");
  }
  Tree->parameter();

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

double FunctionTreeEstimator::evaluate() noexcept {
  auto val = std::dynamic_pointer_cast<Value<double>>(Tree->parameter());
  return val->value();
}

void FunctionTreeEstimator::updateParametersFrom(
    const std::vector<double> &params) {
  size_t pos = 0;
  for (auto p : Parameters.doubleParameters()) {
    p->setValue(params[pos]);
    ++pos;
  }
}

std::vector<ComPWA::Parameter> FunctionTreeEstimator::getParameters() const {
  std::vector<ComPWA::Parameter> params;
  for (auto p : Parameters.doubleParameters()) {
    params.push_back(ComPWA::Parameter{p->name(), p->value()});
  }
  return params;
}

std::string FunctionTreeEstimator::print(int level) const {
  return Tree->print(level);
}

std::shared_ptr<TreeNode> FunctionTreeEstimator::getFunctionTree() const {
  return Tree;
}
ParameterList FunctionTreeEstimator::getParameterList() const {
  return Parameters;
}

FitParameterList
createFitParameterList(ComPWA::FunctionTree::ParameterList Parameters) {
  FitParameterList Pars;
  for (auto x : Parameters.doubleParameters()) {
    ComPWA::FitParameter<double> p;
    p.Value = x->value();
    p.Name = x->name();
    p.HasBounds = x->hasBounds();
    if (p.HasBounds) {
      p.Bounds = x->bounds();
    }
    if (x->hasError()) {
      p.Error = x->error();
    }
    p.IsFixed = x->isFixed();
    Pars.push_back(p);
  }
  return Pars;
}

} // namespace FunctionTree
} // namespace ComPWA
