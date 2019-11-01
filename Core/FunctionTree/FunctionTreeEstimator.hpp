#ifndef CORE_FUNCTIONTREEESTIMATOR_HPP_
#define CORE_FUNCTIONTREEESTIMATOR_HPP_

#include <memory>

#include "Core/FitParameter.hpp"
#include "Core/FunctionTree/ParameterList.hpp"
#include "Estimator/Estimator.hpp"

namespace ComPWA {
namespace FunctionTree {

class FunctionTree;

class FunctionTreeEstimator : public ComPWA::Estimator::Estimator<double> {
public:
  FunctionTreeEstimator(std::shared_ptr<FunctionTree> tree,
                        ParameterList parameters);

  FunctionTreeEstimator(const FunctionTreeEstimator &other) = delete;

  FunctionTreeEstimator(FunctionTreeEstimator &&other) = default;

  double evaluate() noexcept;

  void updateParametersFrom(const std::vector<double> &params);

  std::vector<ComPWA::Parameter> getParameters() const;

  std::string print(int level) const;

  std::shared_ptr<FunctionTree> getFunctionTree() const;
  ParameterList getParameterList() const;

private:
  std::shared_ptr<FunctionTree> Tree;
  ParameterList Parameters;
};

FitParameterList
createFitParameterList(ComPWA::FunctionTree::ParameterList Parameters);

} // namespace FunctionTree
} // namespace ComPWA

#endif
