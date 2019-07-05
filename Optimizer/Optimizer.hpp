// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_OPTIMIZER_HPP_
#define COMPWA_OPTIMIZER_HPP_

#include "Core/FitParameter.hpp"
#include "Estimator/Estimator.hpp"

namespace ComPWA {
namespace Optimizer {

///
/// This class template provides the interface to optimization libraries.
///
/// Note: The dynamic polymorphism by inheriting from Optimizer is not useful
/// within ComPWA, since Optimizers are not never passed to another part of the
/// code.
///
template <typename FitResultType, typename EstimatorType = double>
class Optimizer {
public:
  virtual ~Optimizer() = default;
  /// Finds the optimal value of the Estimator, by varying its parameters
  virtual FitResultType optimize(Estimator::Estimator<EstimatorType> &Estimator,
                                 FitParameterList FitParameters) = 0;
};

} // namespace Optimizer
} // namespace ComPWA

#endif
