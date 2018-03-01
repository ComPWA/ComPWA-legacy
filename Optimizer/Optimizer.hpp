// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Optimizer base class.
///

#ifndef _OPTIMIZER_HPP_
#define _OPTIMIZER_HPP_

#include <vector>
#include <memory>

#include "Core/ParameterList.hpp"
#include "Core/FitResult.hpp"

namespace ComPWA {
namespace Optimizer {

///
/// \class Optimizer
/// This class provides the interface to (external) optimization libraries or
/// routines. As it is pure virtual, one needs at least one implementation to
/// provide an optimizer for the analysis which varies free model-parameters. If
/// a new optimizer is derived from and fulfills this base-class, no change in
/// other modules are necessary to work with the new optimizer library or
/// routine.
///
class Optimizer {

public:
  Optimizer() {}

  virtual ~Optimizer() {  }

  virtual std::shared_ptr<FitResult> exec(ParameterList &par) = 0;
};

} // ns::Optimizer
} // ns::ComPWA

#endif
