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
    : Name(name), UnnormalizedIntensity(intensity), PreviousFitParameters(),
      Integrator(integrator) {

  Normalization = 1.0 / Integrator->integrate(UnnormalizedIntensity);
  UnnormalizedIntensity->addFitParametersTo(PreviousFitParameters);
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
void NormalizationIntensityDecorator::addFitParametersTo(
    std::vector<double> &FitParameters) {
  UnnormalizedIntensity->addFitParametersTo(FitParameters);
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
  std::vector<double> TempParList;
  UnnormalizedIntensity->addFitParametersTo(TempParList);
  for (unsigned int i = 0; i < TempParList.size(); ++i) {
    if (TempParList[i] != PreviousFitParameters[i]) {
      *const_cast<std::vector<double> *>(&PreviousFitParameters) = TempParList;
      return true;
    }
  }
  return false;
}

} // namespace Physics
} // namespace ComPWA
