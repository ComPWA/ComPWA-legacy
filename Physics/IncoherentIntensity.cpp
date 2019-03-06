// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/IncoherentIntensity.hpp"

namespace ComPWA {
namespace Physics {

IncoherentIntensity::IncoherentIntensity(
    const std::string &name,
    const std::vector<std::shared_ptr<ComPWA::Intensity>> &intensities)
    : Name(name), Intensities(intensities) {}

double IncoherentIntensity::evaluate(const ComPWA::DataPoint &point) const {
  double result(0.0);
  for (auto const Intensity : Intensities) {
    result += Intensity->evaluate(point);
  }
  return result;
}

std::shared_ptr<ComPWA::FunctionTree>
IncoherentIntensity::createFunctionTree(const ParameterList &DataSample,
                                        const std::string &suffix) const {

  size_t n = DataSample.mDoubleValue(0)->values().size();

  auto NodeName = "IncoherentIntensity(" + Name + ")" + suffix;

  auto tr = std::make_shared<FunctionTree>(
      NodeName, MDouble("", n), std::make_shared<AddAll>(ParType::MDOUBLE));

  for (auto i : Intensities) {
    auto intensTree = i->createFunctionTree(DataSample, suffix);
    tr->insertTree(intensTree, NodeName);
  }
  return tr;
}

void IncoherentIntensity::addUniqueParametersTo(ParameterList &list) {
  for (auto i : Intensities)
    i->addUniqueParametersTo(list);
}

void IncoherentIntensity::addFitParametersTo(
    std::vector<double> &FitParameters) {
  for (auto i : Intensities) {
    i->addFitParametersTo(FitParameters);
  }
}

void IncoherentIntensity::updateParametersFrom(const ParameterList &list) {
  for (auto i : Intensities)
    i->updateParametersFrom(list);
}

std::vector<std::shared_ptr<ComPWA::Intensity>>
IncoherentIntensity::getIntensities() const {
  return Intensities;
}

} // namespace Physics
} // namespace ComPWA
