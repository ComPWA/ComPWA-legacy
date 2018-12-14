// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_OPTIMIZER_HPP_
#define COMPWA_OPTIMIZER_HPP_

#include <memory>

#include "Core/FitResult.hpp"
#include "Core/ParameterList.hpp"

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

  virtual ~Optimizer() {}

  virtual std::shared_ptr<FitResult> exec(ParameterList &par) = 0;
};

} // namespace Optimizer
} // namespace ComPWA

#endif
