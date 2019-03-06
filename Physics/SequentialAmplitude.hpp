// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_SEQUENTIALAMPLITUDE_HPP_
#define COMPWA_PHYSICS_SEQUENTIALAMPLITUDE_HPP_

#include "Amplitude.hpp"

namespace ComPWA {
namespace Physics {

class SequentialAmplitude : public NamedAmplitude {

public:
  SequentialAmplitude(
      const std::string &name,
      const std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>>
          &PartialAmplitudes_,
      std::complex<double> PreFactor_ = std::complex<double>(1., 0.));

  std::complex<double> evaluate(const DataPoint &point) const final;

  void updateParametersFrom(const ParameterList &list) final;
  void addUniqueParametersTo(ParameterList &list) final;
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<FunctionTree>
  createFunctionTree(const ParameterList &DataSample,
                     const std::string &suffix) const final;

private:
  std::vector<std::shared_ptr<ComPWA::Physics::Amplitude>> PartialAmplitudes;

  /// The PreFactor is a technical feature, which can be used to give an
  /// amplitude an additional factor. This is useful for example in the case
  /// of parity conservation in the helicity formalism, when two amplitudes are
  /// equally strong, but have a pre-factor of -1 with respect to each other.
  std::complex<double> PreFactor;
};

} // namespace Physics
} // namespace ComPWA

#endif
