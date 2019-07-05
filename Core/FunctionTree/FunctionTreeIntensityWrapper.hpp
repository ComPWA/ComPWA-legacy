#ifndef CORE_FUNCTIONTREEINTENSITYWRAPPER_HPP_
#define CORE_FUNCTIONTREEINTENSITYWRAPPER_HPP_

#include <memory>

#include "Core/Function.hpp"
#include "Core/FunctionTree/ParameterList.hpp"

namespace ComPWA {
class Kinematics;
namespace FunctionTree {

class OldIntensity;
class FunctionTree;

class FunctionTreeIntensityWrapper : public Intensity {
public:
  /*FunctionTreeIntensityWrapper(std::shared_ptr<FunctionTree> tree,
                               ParameterList parameters, ParameterList data);*/
  FunctionTreeIntensityWrapper(
      std::shared_ptr<ComPWA::FunctionTree::OldIntensity> oldintens,
      size_t VariableCount, std::string name = "intensity");
  FunctionTreeIntensityWrapper(
      std::shared_ptr<ComPWA::FunctionTree::OldIntensity> oldintens,
      std::shared_ptr<ComPWA::Kinematics> kin, std::string name = "intensity");
  std::vector<double> evaluate(const std::vector<std::vector<double>> &data);
  void updateParametersFrom(const std::vector<double> &params);
  std::vector<double> getParameters() const;

  std::shared_ptr<OldIntensity> getOldIntensity() { return OldIntens; }
  ParameterList getParameterList() { return Parameters; }
  const ParameterList &getUserParameterList() { return UserParameters; }

private:
  void updateDataContainers(const std::vector<std::vector<double>> &data);
  std::shared_ptr<OldIntensity> OldIntens;
  std::shared_ptr<FunctionTree> Tree;
  ParameterList Parameters;
  ParameterList UserParameters;
  ParameterList Data;
};

} // namespace FunctionTree
} // namespace ComPWA

#endif
