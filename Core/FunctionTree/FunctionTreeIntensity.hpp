#ifndef CORE_FUNCTIONTREEINTENSITY_HPP_
#define CORE_FUNCTIONTREEINTENSITY_HPP_

#include <memory>

#include "Core/Function.hpp"
#include "Core/FunctionTree/ParameterList.hpp"
#include "FunctionTreeEstimator.hpp"

namespace ComPWA {
class Kinematics;
namespace FunctionTree {
class TreeNode;

class FunctionTreeIntensity : public Intensity {
public:
  FunctionTreeIntensity(std::shared_ptr<TreeNode> Tree_,
                        ParameterList Parameters_, ParameterList Data_);

  FunctionTreeIntensity(const FunctionTreeIntensity &other) = delete;

  FunctionTreeIntensity(FunctionTreeIntensity &&other) = default;

  std::vector<double>
  evaluate(const std::vector<std::vector<double>> &data) noexcept;

  void updateParametersFrom(const std::vector<double> &params);
  std::vector<ComPWA::Parameter> getParameters() const;

  std::tuple<std::shared_ptr<ComPWA::FunctionTree::TreeNode>,
             ComPWA::FunctionTree::ParameterList>
  bind(const std::vector<std::vector<double>> &data);

  std::string print(int level) const;

private:
  void updateDataContainers(const std::vector<std::vector<double>> &data);

  std::shared_ptr<TreeNode> Tree;
  ParameterList Parameters;
  ParameterList Data;
};

void updateDataContainers(ParameterList Data,
                          const std::vector<std::vector<double>> &data);

} // namespace FunctionTree
} // namespace ComPWA

#endif
