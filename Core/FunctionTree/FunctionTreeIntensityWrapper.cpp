#include "Core/FunctionTree/FunctionTreeIntensityWrapper.hpp"
#include "Core/FunctionTree/FunctionTree.hpp"
#include "Core/FunctionTree/Intensity.hpp"
#include "Core/FunctionTree/Value.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace FunctionTree {
/*FunctionTreeIntensityWrapper::FunctionTreeIntensityWrapper(
    std::shared_ptr<FunctionTree> tree, ParameterList parameters,
    ParameterList data)
    : Tree(tree), Parameters(parameters), Data(data) {}*/

FunctionTreeIntensityWrapper::FunctionTreeIntensityWrapper(
    std::shared_ptr<OldIntensity> oldintens, size_t VariableCount,
    std::string name)
    : OldIntens(oldintens) {
  for (size_t i = 0; i < VariableCount; ++i) {
    std::vector<double> temp;
    Data.addValue(std::make_shared<Value<std::vector<double>>>(temp));
  }
  Tree = oldintens->createFunctionTree(Data, name),
  oldintens->addUniqueParametersTo(Parameters);
  UserParameters.DeepCopy(Parameters);
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

FunctionTreeIntensityWrapper::FunctionTreeIntensityWrapper(
    std::shared_ptr<OldIntensity> oldintens,
    std::shared_ptr<ComPWA::Kinematics> kin, std::string name)
    : FunctionTreeIntensityWrapper(
          oldintens, kin->getKinematicVariableNames().size(), name) {}

std::vector<double> FunctionTreeIntensityWrapper::evaluate(
    const std::vector<std::vector<double>> &data) {
  updateDataContainers(data);
  auto val =
      std::dynamic_pointer_cast<Value<std::vector<double>>>(Tree->parameter());
  return val->value();
}

void FunctionTreeIntensityWrapper::updateDataContainers(
    const std::vector<std::vector<double>> &data) {
  // just loop over the vectors and fill in the data
  if (Data.mDoubleValues().size() > data.size()) {
    std::stringstream ss;
    ss << "FunctionTreeIntensityWrapper::updateDataContainers(): given data "
          "container does not have enough variables! (required: "
       << Data.mDoubleValues().size() << ", given: " << data.size() << ")";
    throw std::out_of_range(ss.str());
  }
  for (size_t i = 0; i < Data.mDoubleValues().size(); ++i) {
    Data.mDoubleValue(i)->setValue(data[i]);
  }
}

void FunctionTreeIntensityWrapper::updateParametersFrom(
    const std::vector<double> &params) {
  size_t pos = 0;
  for (auto p : Parameters.doubleParameters()) {
    p->setValue(params[pos]);
    ++pos;
  }
}

std::vector<double> FunctionTreeIntensityWrapper::getParameters() const {
  std::vector<double> params;
  for (auto p : Parameters.doubleParameters()) {
    params.push_back(p->value());
  }
  return params;
}

} // namespace FunctionTree
} // namespace ComPWA
