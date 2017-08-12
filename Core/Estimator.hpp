// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Estimator Interface Base-Class.
/*! \class Estimator
 * @file Estimator.hpp
 * This class provides the interface to classes which estimate the "closeness"
 * of
 * the modeled intensity to the data. As it is pure virtual, one needs at least
 * one implementation to provide an estimator for the analysis. If a new
 * estimator
 * is derived from and fulfills this base-class, no change in other modules are
 * necessary to work with the new estimation function. As it is derived from
 * ControlParameter, it can be used in the optimizer modules.
*/

#ifndef _ESTIMATOR_HPP_
#define _ESTIMATOR_HPP_

#include <vector>
#include <memory>

#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/AmpIntensity.hpp"

namespace ComPWA {

class IEstimator {

public:
  virtual double controlParameter(ParameterList &minPar) = 0;

  virtual bool HasTree() = 0;

  virtual std::shared_ptr<FunctionTree> GetTree() = 0;

  virtual std::shared_ptr<AmpIntensity> GetIntensity() = 0;

  virtual int nCalls() { return _calls; }

  virtual int getNEvents() { return -999; }

protected:
  std::shared_ptr<FunctionTree> f;
  
  int _calls;
};

} /* namespace ComPWA */

#endif
