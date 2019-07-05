// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_NORMALIZATIONINTENSITYDECORATOR_HPP_
#define PHYSICS_NORMALIZATIONINTENSITYDECORATOR_HPP_

#include "Core/FunctionTree/Intensity.hpp"

namespace ComPWA {

namespace Tools {
class IntegrationStrategy;
}

namespace Physics {

class NormalizationIntensityDecorator
    : public ComPWA::FunctionTree::OldIntensity {
public:
  NormalizationIntensityDecorator(
      const std::string &name,
      std::shared_ptr<ComPWA::FunctionTree::OldIntensity> intensity,
      std::shared_ptr<ComPWA::Tools::IntegrationStrategy> integrator);

  double evaluate(const ComPWA::DataPoint &point) const final;

  void
  updateParametersFrom(const ComPWA::FunctionTree::ParameterList &list) final;
  void addUniqueParametersTo(ComPWA::FunctionTree::ParameterList &list) final;
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createFunctionTree(const ComPWA::FunctionTree::ParameterList &DataSample,
                     const std::string &suffix) const final;

  std::shared_ptr<ComPWA::FunctionTree::OldIntensity>
  getUnnormalizedIntensity() const;

private:
  bool checkParametersChanged() const;

  std::string Name;
  std::shared_ptr<ComPWA::FunctionTree::OldIntensity> UnnormalizedIntensity;

  double Normalization;
  std::vector<double> PreviousFitParameters;
  /// Phsp sample for numerical integration
  std::shared_ptr<ComPWA::Tools::IntegrationStrategy> Integrator;
};

} // namespace Physics
} // namespace ComPWA

#endif
