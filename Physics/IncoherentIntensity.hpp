// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_INCOHERENT_INTENSITY_HPP_
#define COMPWA_PHYSICS_INCOHERENT_INTENSITY_HPP_

#include "Core/Intensity.hpp"
//#include "Core/Properties.hpp"

namespace ComPWA {

class Kinematics;

namespace Physics {

class IncoherentIntensity
    : public ComPWA::Intensity,
      public std::enable_shared_from_this<IncoherentIntensity> {

public:
  IncoherentIntensity(
      const std::string &name,
      const std::vector<std::shared_ptr<ComPWA::Intensity>> &intensities);

  double evaluate(const ComPWA::DataPoint &point) const final;

  void updateParametersFrom(const ParameterList &list) final;
  void addUniqueParametersTo(ParameterList &list) final;
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<FunctionTree>
  createFunctionTree(const ParameterList &DataSample,
                     const std::string &suffix) const final;

  std::vector<std::shared_ptr<ComPWA::Intensity>> getIntensities() const;

private:
  std::string Name;
  std::vector<std::shared_ptr<ComPWA::Intensity>> Intensities;
};

} // namespace Physics
} // namespace ComPWA

#endif
