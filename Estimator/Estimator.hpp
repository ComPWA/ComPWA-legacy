// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_ESTIMATOR_ESTIMATOR_HPP_
#define COMPWA_ESTIMATOR_ESTIMATOR_HPP_

namespace ComPWA {
class ParameterList;
namespace Estimator {

///
/// This class provides the interface to classes which estimate the "closeness"
/// of the modeled intensity to the data. Any derived Estimator can be used with
/// any derived Optimizer.
///
class Estimator {

public:
  virtual ~Estimator(){};

  /// Evaluates the Estimator, which calculates the "distance" of the
  /// Intensity from the DataPoints (or more generally a model from the data).
  /// The Optimizer tries to minimize/optimize the returned value of the
  /// Estimators evaluate function.
  virtual double evaluate() const = 0;
};

} // namespace Estimator
} // namespace ComPWA

#endif
