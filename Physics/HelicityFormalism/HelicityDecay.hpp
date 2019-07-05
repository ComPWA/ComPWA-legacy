// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_HELICITYFORMALISM_HELICITYDECAY_HPP_
#define COMPWA_PHYSICS_HELICITYFORMALISM_HELICITYDECAY_HPP_

#include <memory>

#include "Physics/Amplitude.hpp"
#include "Physics/Dynamics/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

///
/// \class HelicityDecay
/// Represents a two-body decay within the helicity formalism.
///
class HelicityDecay : public NamedAmplitude {

public:
  HelicityDecay(
      const std::string &name,
      std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD>
          AngularFunction_,
      std::shared_ptr<ComPWA::Physics::Dynamics::AbstractDynamicalFunction>
          DynamicFunction_,
      unsigned int datapos, double prefactor = 1.0);

  std::complex<double> evaluate(const DataPoint &point) const final;

  void
  updateParametersFrom(const ComPWA::FunctionTree::ParameterList &list) final;
  void addUniqueParametersTo(ComPWA::FunctionTree::ParameterList &list) final;
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createFunctionTree(const ComPWA::FunctionTree::ParameterList &DataSample,
                     const std::string &suffix) const final;

private:
  // ComPWA::Physics::SubSystem SubSys;

  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD>
      AngularFunction;

  std::shared_ptr<ComPWA::Physics::Dynamics::AbstractDynamicalFunction>
      DynamicFunction;

  /// Position where variables are stored in DataPoint.
  /// We expect to find the invariant mass of the system at @param DataPosition,
  /// theta at @param DataPosition+1 and phi at @param DataPosition+2
  unsigned int DataPosition;

  double PreFactor;
};

} // namespace HelicityFormalism
} // namespace Physics
} // namespace ComPWA

#endif
