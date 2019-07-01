#ifndef CORE_FUNCTIONTREEINTENSITYWRAPPER_HPP_
#define CORE_FUNCTIONTREEINTENSITYWRAPPER_HPP_

#include <memory>

#include "Core/Function.hpp"
#include "Core/ParameterList.hpp"

namespace ComPWA {

class Kinematics;
class OldIntensity;
class FunctionTree;

class FunctionTreeIntensityWrapper : public Intensity {
public:
  /*FunctionTreeIntensityWrapper(std::shared_ptr<FunctionTree> tree,
                               ParameterList parameters, ParameterList data);*/
  FunctionTreeIntensityWrapper(std::shared_ptr<ComPWA::OldIntensity> oldintens,
                               size_t VariableCount,
                               std::string name = "intensity");
  FunctionTreeIntensityWrapper(std::shared_ptr<ComPWA::OldIntensity> oldintens,
                               std::shared_ptr<ComPWA::Kinematics> kin,
                               std::string name = "intensity");
  double evaluate(const std::vector<double> &data);
  void updateParametersFrom(const std::vector<double> &params);
  std::vector<double> getParameters() const;

private:
  void updateDataContainers(const std::vector<double> &data);
  std::shared_ptr<FunctionTree> Tree;
  ParameterList Parameters;
  ParameterList Data;
};

} /* namespace ComPWA */

#endif /* CORE_FUNCTIONTREEINTENSITYWRAPPER_HPP_ */
