// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef OPTIMIZER_MINUIT2_MINUITIF_HPP
#define OPTIMIZER_MINUIT2_MINUITIF_HPP

#include "Core/FitParameter.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"
#include "Optimizer/Optimizer.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

///
/// \class MinuitIF
/// Wrapper of the Minuit2 Optimizer library. This class provides a wrapper
/// around the Minuit2 library. It fulfills the
/// Optimizer interface to be easily adapted to other modules. The data needs to
/// be provided with the ControlParameter interface.
///
class MinuitIF : public Optimizer<MinuitResult> {
public:
  MinuitIF() { setStrategy("medium"); }
  MinuitResult optimize(ComPWA::Estimator::Estimator<double> &Estimator,
                        ComPWA::FitParameterList InitialParameters);

  bool UseHesse = 1;
  bool UseMinos = 0;

  /// Minuit strategy (low, medium(default), high)
  /// See
  /// https://root.cern.ch/root/htmldoc/guides/minuit2/Minuit2.html#m-strategy
  /// Sets Minuit configuration variables to pre-defined values
  void setStrategy(std::string strategy);

  // Minuit2 configuration variables
  unsigned int GradientNCycles;
  double GradientStepTolerance;
  double GradientTolerance;
  unsigned int HessianNCycles;
  unsigned int HessianGradientNCycles;
  double HessianStepTolerance;
  double HessianG2Tolerance;

private:
  /// Check if a pre-defined strategy is set or if custom settings are used
  std::string checkStrategy();
};

} // namespace Minuit2
} // namespace Optimizer
} // namespace ComPWA

#endif
