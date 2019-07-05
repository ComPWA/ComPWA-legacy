// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>

#include "Core/Event.hpp"
#include "NormalizationAmplitudeDecorator.hpp"
#include "Physics/CoherentIntensity.hpp"
#include "Tools/Integration.hpp"

namespace ComPWA {
namespace Physics {

using namespace ComPWA::FunctionTree;

NormalizationAmplitudeDecorator::NormalizationAmplitudeDecorator(
    const std::string &name, std::shared_ptr<NamedAmplitude> amplitude,
    std::shared_ptr<ComPWA::Tools::IntegrationStrategy> integrator)
    : NamedAmplitude(name), UnnormalizedAmplitude(amplitude),
      NormedAmplitude(std::make_shared<CoherentIntensity>(
          amplitude->getName(),
          std::vector<std::shared_ptr<NamedAmplitude>>{amplitude})),
      PreviousFitParameters(), Integrator(integrator) {
  // Normalization = std::sqrt(1.0 / Integrator->integrate(NormedAmplitude));
  // UnnormalizedAmplitude->addFitParametersTo(PreviousFitParameters);
}

std::complex<double> NormalizationAmplitudeDecorator::evaluate(
    const ComPWA::DataPoint &point) const {

  /* double Norm(Normalization);
   if (checkParametersChanged()) {
     LOG(DEBUG) << "NormalizationAmplitudeDecorator::evaluate(): recalculating "
                   "normalization for amplitude";
     Norm = 1.0 / Integrator->integrate(NormedAmplitude);
   }*/
  double Norm(1.0);

  return Norm * UnnormalizedAmplitude->evaluate(point);
}

void NormalizationAmplitudeDecorator::updateParametersFrom(
    const ParameterList &list) {
  UnnormalizedAmplitude->updateParametersFrom(list);
}

void NormalizationAmplitudeDecorator::addUniqueParametersTo(
    ParameterList &list) {
  UnnormalizedAmplitude->addUniqueParametersTo(list);
}
void NormalizationAmplitudeDecorator::addFitParametersTo(
    std::vector<double> &FitParameters) {
  UnnormalizedAmplitude->addFitParametersTo(FitParameters);
}

std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
NormalizationAmplitudeDecorator::createFunctionTree(
    const ParameterList &DataSample, const std::string &suffix) const {

  auto NodeName = "NormalizedAmplitude(" + getName() + ")" + suffix;
  size_t n = DataSample.mDoubleValue(0)->values().size();

  auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
      NodeName, MComplex("", n), std::make_shared<MultAll>(ParType::MCOMPLEX));

  tr->createNode("sqrt(Normalization)", ValueFactory(ParType::DOUBLE),
                 std::make_shared<SquareRoot>(ParType::DOUBLE), NodeName);
  auto normtree = Integrator->createFunctionTree(NormedAmplitude, "_ampnorm");
  tr->insertTree(normtree, "sqrt(Normalization)");
  auto intenstree =
      UnnormalizedAmplitude->createFunctionTree(DataSample, suffix);
  tr->insertTree(intenstree, NodeName);

  return tr;
}

std::shared_ptr<const Amplitude>
NormalizationAmplitudeDecorator::getUnnormalizedAmplitude() const {
  return UnnormalizedAmplitude;
}

bool NormalizationAmplitudeDecorator::checkParametersChanged() const {
  std::vector<double> TempParList;
  UnnormalizedAmplitude->addFitParametersTo(TempParList);
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
