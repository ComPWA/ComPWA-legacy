#ifndef CORE_FUNCTIONTREEESTIMATORWRAPPER_HPP_
#define CORE_FUNCTIONTREEESTIMATORWRAPPER_HPP_

#include <memory>

#include "Core/FunctionTree/ParameterList.hpp"
#include "Estimator/Estimator.hpp"

namespace ComPWA {
namespace FunctionTree {

class FunctionTree;

class FunctionTreeEstimatorWrapper
    : public ComPWA::Estimator::Estimator<double> {
public:
  FunctionTreeEstimatorWrapper(std::shared_ptr<FunctionTree> tree,
                               ParameterList parameters,
                               ParameterList userparameters);
  double evaluate();
  void updateParametersFrom(const std::vector<double> &params);
  std::vector<double> getParameters() const;

  std::shared_ptr<FunctionTree> getFunctionTree() { return Tree; }
  ParameterList &getParameterList() { return Parameters; };
  ParameterList &getUserParameterList() { return UserParameters; };

  std::string print(int level) const;

private:
  std::shared_ptr<FunctionTree> Tree;
  ParameterList Parameters;
  ParameterList UserParameters;
};

} // namespace FunctionTree
} /* namespace ComPWA */

#endif
