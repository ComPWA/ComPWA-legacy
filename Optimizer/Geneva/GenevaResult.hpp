// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_OPTIMIZER_GENEVARESULT_HPP_
#define COMPWA_OPTIMIZER_GENEVARESULT_HPP_

#include "Core/FitResult.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Geneva {

struct GenevaResult : public FitResult {
  GenevaResult() = default;
  GenevaResult(FitResult Result) : FitResult(Result) {}
  // TODO: add more properties here
};

} // namespace Geneva
} // namespace Optimizer
} // namespace ComPWA

#endif
