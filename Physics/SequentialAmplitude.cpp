// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "SequentialAmplitude.hpp"

namespace ComPWA {
namespace Physics {

SequentialAmplitude::SequentialAmplitude(
    const std::string &name,
    const std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>>
        &PartialAmplitudes_,
    std::complex<double> PreFactor_)
    : NamedAmplitude(name), PartialAmplitudes(PartialAmplitudes_),
      PreFactor(PreFactor_) {}

std::complex<double>
SequentialAmplitude::evaluate(const DataPoint &point) const {
  std::complex<double> result(PreFactor);
  for (auto i : PartialAmplitudes)
    result *= i->evaluate(point);

  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return result;
};

std::shared_ptr<ComPWA::FunctionTree>
SequentialAmplitude::createFunctionTree(const ParameterList &DataSample,
                                        const std::string &suffix) const {

  size_t n = DataSample.mDoubleValue(0)->values().size();

  auto NodeName = "SequentialPartialAmplitude(" + getName() + ")" + suffix;
  auto tr = std::make_shared<FunctionTree>(
      NodeName, MComplex("", n), std::make_shared<MultAll>(ParType::MCOMPLEX));
  tr->createLeaf("Prefactor", PreFactor, NodeName);

  for (auto i : PartialAmplitudes) {
    std::shared_ptr<FunctionTree> resTree =
        i->createFunctionTree(DataSample, suffix);
    if (!resTree->sanityCheck())
      throw std::runtime_error("SequentialAmplitude::createFunctionTree : tree "
                               "didn't pass sanity check!");
    resTree->parameter();
    tr->insertTree(resTree, NodeName);
  }

  return tr;
}

void SequentialAmplitude::addUniqueParametersTo(ParameterList &list) {
  for (auto i : PartialAmplitudes)
    i->addUniqueParametersTo(list);
}
void SequentialAmplitude::addFitParametersTo(
    std::vector<double> &FitParameters) {
  for (auto i : PartialAmplitudes)
    i->addFitParametersTo(FitParameters);
}

void SequentialAmplitude::updateParametersFrom(const ParameterList &list) {
  for (auto i : PartialAmplitudes)
    i->updateParametersFrom(list);
}

} // namespace Physics
} // namespace ComPWA
