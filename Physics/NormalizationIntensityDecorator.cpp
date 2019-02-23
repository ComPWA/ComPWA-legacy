// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "NormalizationIntensityDecorator.hpp"
#include "Core/Event.hpp"
#include "Tools/Integration.hpp"

namespace ComPWA {
namespace Physics {

NormalizationIntensityDecorator::NormalizationIntensityDecorator(
    const std::string &name, std::shared_ptr<ComPWA::Intensity> intensity,
    std::shared_ptr<ComPWA::Tools::IntegrationStrategy> integrator)
    : Name(name), UnnormalizedIntensity(intensity), Integrator(integrator) {

  Normalization = 1.0 / Integrator->integrate(UnnormalizedIntensity);
  ParameterList TempParList;
  UnnormalizedIntensity->addUniqueParametersTo(TempParList);
  PreviousParameterList.DeepCopy(TempParList);
}

double NormalizationIntensityDecorator::evaluate(
    const ComPWA::DataPoint &point) const {

  double Norm(Normalization);
  if (checkParametersChanged()) {
    LOG(DEBUG) << "NormalizationIntensityDecorator::evaluate(): recalculating "
                  "normalization for intensity";
    Norm = 1.0 / Integrator->integrate(UnnormalizedIntensity);
  }

  return Norm * UnnormalizedIntensity->evaluate(point);
}

void NormalizationIntensityDecorator::updateParametersFrom(
    const ParameterList &list) {
  UnnormalizedIntensity->updateParametersFrom(list);
}

void NormalizationIntensityDecorator::addUniqueParametersTo(
    ParameterList &list) {
  UnnormalizedIntensity->addUniqueParametersTo(list);
}

std::shared_ptr<FunctionTree>
NormalizationIntensityDecorator::createFunctionTree(
    const ParameterList &DataSample, const std::string &suffix) const {

  auto NodeName = "NormalizedIntensity(" + Name + ")" + suffix;
  size_t n = DataSample.mDoubleValue(0)->values().size();

  auto tr = std::make_shared<FunctionTree>(
      NodeName, MDouble("", n), std::make_shared<MultAll>(ParType::MDOUBLE));

  auto normtree =
      Integrator->createFunctionTree(UnnormalizedIntensity, "_norm");
  tr->insertTree(normtree, NodeName);
  auto intenstree =
      UnnormalizedIntensity->createFunctionTree(DataSample, suffix);
  tr->insertTree(intenstree, NodeName);

  return tr;
}

std::shared_ptr<const ComPWA::Intensity>
NormalizationIntensityDecorator::getUnnormalizedIntensity() const {
  return UnnormalizedIntensity;
}

bool NormalizationIntensityDecorator::checkParametersChanged() const {
  ParameterList TempParList;
  UnnormalizedIntensity->addUniqueParametersTo(TempParList);
  for (unsigned int i = 0; i < TempParList.doubleParameters().size(); ++i) {
    if (TempParList.doubleParameter(i)->value() !=
        PreviousParameterList.doubleParameter(i)->value()) {
      const_cast<ParameterList *>(&PreviousParameterList)
          ->DeepCopy(TempParList);
      return true;
    }
  }
  return false;
}

} // namespace Physics
} // namespace ComPWA
