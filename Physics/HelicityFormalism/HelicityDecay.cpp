// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/HelicityFormalism/HelicityDecay.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

HelicityDecay::HelicityDecay(
    const std::string &name,
    std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD>
        AngularFunction_,
    std::shared_ptr<ComPWA::Physics::Dynamics::AbstractDynamicalFunction>
        DynamicFunction_,
    unsigned int datapos, double prefactor)
    : NamedAmplitude(name), AngularFunction(AngularFunction_),
      DynamicFunction(DynamicFunction_), DataPosition(datapos),
      PreFactor(prefactor) {}

std::complex<double> HelicityDecay::evaluate(const DataPoint &point) const {
  std::complex<double> result(PreFactor, 0);
  result *=
      AngularFunction->evaluate(point, DataPosition + 1, DataPosition + 2);
  result *= DynamicFunction->evaluate(point, DataPosition);

  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return result;
};

std::shared_ptr<FunctionTree>
HelicityDecay::createFunctionTree(const ParameterList &DataSample,
                                  const std::string &suffix) const {

  size_t n = DataSample.mDoubleValue(0)->values().size();

  std::string nodeName = "PartialAmplitude(" + getName() + ")" + suffix;

  auto tr = std::make_shared<FunctionTree>(
      nodeName, MComplex("", n), std::make_shared<MultAll>(ParType::MCOMPLEX));
  tr->createLeaf("PreFactor", PreFactor, nodeName);
  tr->insertTree(
      AngularFunction->tree(DataSample, DataPosition + 1, DataPosition + 2),
      nodeName);
  tr->insertTree(
      DynamicFunction->createFunctionTree(DataSample, DataPosition, ""),
      nodeName);

  tr->parameter();
  return tr;
}

void HelicityDecay::addUniqueParametersTo(ParameterList &list) {
  DynamicFunction->addUniqueParametersTo(list);
}
void HelicityDecay::addFitParametersTo(std::vector<double> &FitParameters) {
  DynamicFunction->addFitParametersTo(FitParameters);
}

void HelicityDecay::updateParametersFrom(const ParameterList &list) {
  DynamicFunction->updateParametersFrom(list);
}

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA
