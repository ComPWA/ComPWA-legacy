#ifndef CORE_FUNCTIONTREEINTENSITY_HPP_
#define CORE_FUNCTIONTREEINTENSITY_HPP_

#include <memory>

#include "Core/Function.hpp"
#include "Core/FunctionTree/ParameterList.hpp"

namespace ComPWA {
namespace Data {
struct DataMap;
}
namespace FunctionTree {
class TreeNode;

class FunctionTreeIntensity : public Intensity {
public:
  FunctionTreeIntensity(std::shared_ptr<TreeNode> Tree_,
                        ParameterList Parameters_, ParameterList Data_);

  FunctionTreeIntensity(const FunctionTreeIntensity &other) = delete;

  FunctionTreeIntensity(FunctionTreeIntensity &&other) = default;

  std::vector<double> evaluate(const ComPWA::DataMap &data) noexcept;

  void updateParametersFrom(const std::vector<double> &params);
  std::vector<ComPWA::Parameter> getParameters() const;

  std::tuple<std::shared_ptr<ComPWA::FunctionTree::TreeNode>,
             ComPWA::FunctionTree::ParameterList>
  bind(const ComPWA::DataMap &data);

  std::string print(int level) const;

private:
  void updateDataContainers(const ComPWA::DataMap &data);

  std::shared_ptr<TreeNode> Tree;
  ParameterList Parameters;
  ParameterList Data;
};

void updateDataContainers(ParameterList Data, const ComPWA::DataMap &data);

} // namespace FunctionTree
} // namespace ComPWA

#endif
