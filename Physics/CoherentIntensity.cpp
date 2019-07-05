// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Physics/CoherentIntensity.hpp"
#include "Physics/Amplitude.hpp"

namespace ComPWA {
namespace Physics {

using namespace ComPWA::FunctionTree;

CoherentIntensity::CoherentIntensity(
    const std::string &name,
    const std::vector<std::shared_ptr<ComPWA::Physics::NamedAmplitude>>
        &amplitudes)
    : Name(name), Amplitudes(amplitudes) {}

double CoherentIntensity::evaluate(const DataPoint &point) const {
  std::complex<double> result(0., 0.);
  for (auto i : Amplitudes)
    result += i->evaluate(point);
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return std::norm(result);
};

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
CoherentIntensity::createFunctionTree(const ParameterList &DataSample,
                                      const std::string &suffix) const {

  size_t n = DataSample.mDoubleValue(0)->values().size();

  auto NodeName = "CoherentIntensity(" + Name + ")" + suffix;

  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MDouble("", n), std::make_shared<AbsSquare>(ParType::MDOUBLE));

  tr->createNode("SumOfAmplitudes", MComplex("", n),
                 std::make_shared<AddAll>(ParType::MCOMPLEX), NodeName);

  for (auto i : Amplitudes) {
    std::shared_ptr<ComPWA::FunctionTree::FunctionTree> resTree =
        i->createFunctionTree(DataSample, suffix);
    if (!resTree->sanityCheck())
      throw std::runtime_error("CoherentIntensity::createFunctionTree(): tree "
                               "didn't pass sanity check!");
    resTree->parameter();
    tr->insertTree(resTree, "SumOfAmplitudes");
  }

  return tr;
}

void CoherentIntensity::addUniqueParametersTo(ParameterList &list) {
  for (auto i : Amplitudes) {
    i->addUniqueParametersTo(list);
  }
}

void CoherentIntensity::addFitParametersTo(std::vector<double> &FitParameters) {
  for (auto i : Amplitudes) {
    i->addFitParametersTo(FitParameters);
  }
}

void CoherentIntensity::updateParametersFrom(const ParameterList &list) {
  for (auto i : Amplitudes)
    i->updateParametersFrom(list);
}

const std::vector<std::shared_ptr<ComPWA::Physics::NamedAmplitude>> &
CoherentIntensity::getAmplitudes() const {
  return Amplitudes;
}

} // namespace Physics
} // namespace ComPWA
