// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "StrengthIntensityDecorator.hpp"

namespace ComPWA {
namespace Physics {

StrengthIntensityDecorator::StrengthIntensityDecorator(
    const std::string &name, std::shared_ptr<ComPWA::OldIntensity> Intensity,
    std::shared_ptr<ComPWA::FitParameter> strength)
    : Name(name), UndecoratedIntensity(Intensity), Strength(strength) {}

double StrengthIntensityDecorator::evaluate(const DataPoint &point) const {
  return Strength->value() * UndecoratedIntensity->evaluate(point);
};

void StrengthIntensityDecorator::addUniqueParametersTo(
    ComPWA::ParameterList &list) {
  Strength = list.addUniqueParameter(Strength);
  UndecoratedIntensity->addUniqueParametersTo(list);
}

void StrengthIntensityDecorator::addFitParametersTo(
    std::vector<double> &FitParameters) {
  FitParameters.push_back(Strength->value());
  UndecoratedIntensity->addFitParametersTo(FitParameters);
}

void StrengthIntensityDecorator::updateParametersFrom(
    const ParameterList &list) {
  std::shared_ptr<FitParameter> p;
  try {
    p = FindParameter(Strength->name(), list);
  } catch (std::exception &ex) {
  }
  if (p)
    Strength->updateParameter(p);
  UndecoratedIntensity->updateParametersFrom(list);
}

std::shared_ptr<FunctionTree> StrengthIntensityDecorator::createFunctionTree(
    const ParameterList &DataSample, const std::string &suffix) const {

  size_t n = DataSample.mDoubleValue(0)->values().size();

  auto NodeName = "StrengthIntensity(" + Name + ")" + suffix;

  auto tr = std::make_shared<FunctionTree>(
      NodeName, MDouble("", n), std::make_shared<MultAll>(ParType::MDOUBLE));

  tr->createLeaf("Strength", Strength, NodeName);

  std::shared_ptr<ComPWA::FunctionTree> x =
      UndecoratedIntensity->createFunctionTree(DataSample, "");

  x->parameter();
  tr->insertTree(x, NodeName);

  if (!tr->sanityCheck())
    throw std::runtime_error(
        "StrengthIntensityDecorator::createFunctionTree() | "
        "Tree didn't pass sanity check!");

  return tr;
}

} // namespace Physics
} // namespace ComPWA
