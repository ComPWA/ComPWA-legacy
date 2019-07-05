// Copyright (c) 2014, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "MinuitResult.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

void MinuitResult::print(std::ostream &os) const {
  // print base infos
  FitResult::print(os);

  os << "--------------MINUIT2 FIT INFOS----------------\n";
  if (!IsValid)
    os << "    *** MINIMUM NOT VALID! ***\n";
  os << "Estimated distance to minimum: " << Edm << "\n";
  if (EdmAboveMax)
    os << "    *** EDM IS ABOVE MAXIMUM! ***\n";
  os << "Error definition: " << ErrorDef << "\n";
  os << "Number of calls: " << NFcn << "\n";
  if (HasReachedCallLimit)
    os << "    *** LIMIT OF MAX CALLS REACHED! ***\n";
  if (!HasValidParameters)
    os << "    *** NO VALID SET OF PARAMETERS! ***\n";
  if (!HasValidCov)
    os << "    *** COVARIANCE MATRIX NOT VALID! ***\n";
  if (!HasAccCov)
    os << "    *** COVARIANCE MATRIX NOT ACCURATE! ***\n";
  if (!CovPosDef)
    os << "    *** COVARIANCE MATRIX NOT POSITIVE DEFINITE! ***\n";
  if (HesseFailed)
    os << "    *** HESSE FAILED! ***\n";
  os << "-----------------------------------------------\n";
}

std::ostream &operator<<(std::ostream &os, const MinuitResult &Result) {
  Result.print(os);
  return os;
}

} // namespace Minuit2
} // namespace Optimizer
} // namespace ComPWA
