// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_STRENGTHINTENSITYDECORATOR_HPP_
#define COMPWA_STRENGTHINTENSITYDECORATOR_HPP_

#include "Core/Intensity.hpp"

namespace ComPWA {
namespace Physics {

class StrengthIntensityDecorator : public Intensity {
public:
  StrengthIntensityDecorator(const std::string &name,
                             std::shared_ptr<ComPWA::Intensity> Intensity,
                             std::shared_ptr<ComPWA::FitParameter> strength);

  virtual ~StrengthIntensityDecorator() = default;

  double evaluate(const ComPWA::DataPoint &point) const final;

  void updateParametersFrom(const ParameterList &list) final;
  void addUniqueParametersTo(ParameterList &list) final;
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<FunctionTree>
  createFunctionTree(const ParameterList &DataSample,
                     const std::string &suffix) const final;

private:
  std::string Name;
  std::shared_ptr<ComPWA::Intensity> UndecoratedIntensity;
  std::shared_ptr<ComPWA::FitParameter> Strength;
};

} // namespace Physics
} // namespace ComPWA

#endif
