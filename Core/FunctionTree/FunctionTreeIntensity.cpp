#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/FunctionTree/TreeNode.hpp"
#include "Core/FunctionTree/Value.hpp"
#include "Data/DataSet.hpp"

namespace ComPWA {
namespace FunctionTree {
FunctionTreeIntensity::FunctionTreeIntensity(std::shared_ptr<TreeNode> Tree_,
                                             ParameterList Parameters_,
                                             ParameterList Data_)
    : Tree(Tree_), Parameters(Parameters_), Data(Data_) {
  Tree->parameter();
}

std::vector<double>
FunctionTreeIntensity::evaluate(const ComPWA::DataMap &data) noexcept {
  updateDataContainers(data);
  auto val =
      std::dynamic_pointer_cast<Value<std::vector<double>>>(Tree->parameter());
  return val->value();
}

void FunctionTreeIntensity::updateDataContainers(const ComPWA::DataMap &data) {
  ComPWA::FunctionTree::updateDataContainers(Data, data);
}

void FunctionTreeIntensity::updateParametersFrom(
    const std::vector<double> &params) {
  size_t pos = 0;
  for (auto p : Parameters.doubleParameters()) {
    p->setValue(params[pos]);
    ++pos;
  }
}

std::vector<ComPWA::Parameter> FunctionTreeIntensity::getParameters() const {
  std::vector<ComPWA::Parameter> params;
  for (auto p : Parameters.doubleParameters()) {
    params.push_back(ComPWA::Parameter{p->name(), p->value()});
  }
  return params;
}

std::tuple<std::shared_ptr<ComPWA::FunctionTree::TreeNode>,
           ComPWA::FunctionTree::ParameterList>
FunctionTreeIntensity::bind(const ComPWA::DataMap &data) {
  updateDataContainers(data);
  return std::make_tuple(Tree, Parameters);
}

std::string FunctionTreeIntensity::print(int level) const {
  return Tree->print(level);
}

void updateDataContainers(ParameterList Data, const ComPWA::DataMap &data) {
  // just loop over the vectors and fill in the data
  if (Data.mDoubleValues().size() > data.size()) {
    std::stringstream ss;
    ss << "FunctionTreeIntensity::updateDataContainers(): given data "
          "container does not have enough variables! (required: "
       << Data.mDoubleValues().size() << ", given: " << data.size() << ")";
    throw std::out_of_range(ss.str());
  }
  for (size_t i = 0; i < Data.mDoubleValues().size(); ++i) {
    Data.mDoubleValue(i)->setValue(data.at(Data.mDoubleValue(i)->name()));
  }
}

} // namespace FunctionTree
} // namespace ComPWA
