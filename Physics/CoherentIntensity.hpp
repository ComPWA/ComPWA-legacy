// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_COHERENTAMPLITUDE_HPP_
#define COMPWA_PHYSICS_COHERENTAMPLITUDE_HPP_

#include "Core/Intensity.hpp"

namespace ComPWA {
namespace Physics {

class Amplitude;

class CoherentIntensity
    : public ComPWA::Intensity,
      public std::enable_shared_from_this<CoherentIntensity> {

public:
  CoherentIntensity(
      const std::string &name,
      const std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>>
          &amplitudes);

  virtual ~CoherentIntensity() = default;

  double evaluate(const ComPWA::DataPoint &point) const final;

  void updateParametersFrom(const ParameterList &list) final;
  void addUniqueParametersTo(ParameterList &list) final;

  std::shared_ptr<FunctionTree>
  createFunctionTree(const ParameterList &DataSample,
                     const std::string &suffix) const final;

  const std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>> &
  getAmplitudes() const;

private:
  std::string Name;
  std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>> Amplitudes;
};

} // namespace Physics
} // namespace ComPWA

#endif
