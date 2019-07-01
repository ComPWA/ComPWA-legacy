#include "FunctionTreeIntensityWrapper.hpp"
#include "Core/Intensity.hpp"
#include "FunctionTree.hpp"
#include "Kinematics.hpp"
#include "Value.hpp"

namespace ComPWA {

/*FunctionTreeIntensityWrapper::FunctionTreeIntensityWrapper(
    std::shared_ptr<FunctionTree> tree, ParameterList parameters,
    ParameterList data)
    : Tree(tree), Parameters(parameters), Data(data) {}*/

FunctionTreeIntensityWrapper::FunctionTreeIntensityWrapper(
    std::shared_ptr<ComPWA::OldIntensity> oldintens, size_t VariableCount,
    std::string name) {
  for (size_t i = 0; i < VariableCount; ++i) {
    std::vector<double> temp = {0.0};
    Data.addValue(std::make_shared<Value<std::vector<double>>>(temp));
  }
  Tree = oldintens->createFunctionTree(Data, name),
  oldintens->addUniqueParametersTo(Parameters);
}

FunctionTreeIntensityWrapper::FunctionTreeIntensityWrapper(
    std::shared_ptr<ComPWA::OldIntensity> oldintens,
    std::shared_ptr<ComPWA::Kinematics> kin, std::string name)
    : FunctionTreeIntensityWrapper(
          oldintens, kin->getKinematicVariableNames().size(), name) {}

double FunctionTreeIntensityWrapper::evaluate(const std::vector<double> &data) {
  updateDataContainers(data);
  Tree->UpdateAll(Tree->head());
  auto val =
      std::dynamic_pointer_cast<Value<std::vector<double>>>(Tree->parameter());
  return val->value()[0];
}

void FunctionTreeIntensityWrapper::updateDataContainers(
    const std::vector<double> &data) {
  // just loop over the vectors and fill in the data
  if (Data.mDoubleValues().size() != data.size()) {
    std::stringstream ss;
    ss << "FunctionTreeIntensityWrapper::updateDataContainers(): data "
          "containers do not match in size! (Data: "
       << Data.mDoubleValues().size() << ", data: " << data.size() << ")";
    throw std::out_of_range(ss.str());
  }
  for (size_t i = 0; i < data.size(); ++i) {
    Data.mDoubleValue(i)->values()[0] = data[i];
  }
}

void FunctionTreeIntensityWrapper::updateParametersFrom(
    const std::vector<double> &params) {
  size_t pos = 0;
  for (auto p : Parameters.doubleParameters()) {
    if (!p->isFixed()) {
      p->setValue(params[pos]);
    }
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

} /* namespace ComPWA */
