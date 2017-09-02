// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Estimator interface
///

#ifndef _ESTIMATOR_HPP_
#define _ESTIMATOR_HPP_

#include <vector>
#include <memory>

#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/AmpIntensity.hpp"

namespace ComPWA {

///
/// \class Estimator
/// Estimator Interface Base-Class.
/// This class provides the interface to classes which estimate the "closeness"
/// of
/// the modeled intensity to the data. As it is pure virtual, one needs at least
/// one implementation to provide an estimator for the analysis. If a new
/// estimator
/// is derived from and fulfills this base-class, no change in other modules are
/// necessary to work with the new estimation function. As it is derived from
/// ControlParameter, it can be used in the optimizer modules.
///
class IEstimator {

public:
  /// Calculate value of Estimator given parameters \p par
  virtual double controlParameter(ParameterList &par) = 0;

  /// Get the FunctionTree.
  /// If no FunctionTree is available an exception is thrown.
  virtual std::shared_ptr<FunctionTree> GetTree() = 0;

  /// Kind for status flag during the minimization process.
  /// (e.g. number of likelihood evaluations)
  virtual int NSteps() = 0 ;

protected:
  std::shared_ptr<FunctionTree> f;
};

} // namespace ComPWA

#endif
