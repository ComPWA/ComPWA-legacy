// Copyright (c) 2014, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "MinuitResult.hpp"

#include "Core/Logging.hpp"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

std::vector<double>
getGlobalCorrelations(const ROOT::Minuit2::MnUserParameterState &minState) {
  std::vector<double> GlobalCC;
  if (minState.HasGlobalCC()) {
    GlobalCC = minState.GlobalCC().GlobalCC();
  } else {
    auto NumFreeParameter = minState.Parameters().Trafo().VariableParameters();
    GlobalCC = std::vector<double>(NumFreeParameter, 0);
    LOG(ERROR)
        << "getGlobalCorrelations() | no valid global correlation available!";
  }
  return GlobalCC;
}

MinuitResult::MinuitResult(const FitResult &Result,
                           const ROOT::Minuit2::FunctionMinimum &FMin)
    : FitResult(Result), CovPosDef(FMin.HasPosDefCovar()),
      HasValidParameters(FMin.HasValidParameters()),
      HasValidCov(FMin.HasValidCovariance()),
      HasAccCov(FMin.HasAccurateCovar()),
      HasReachedCallLimit(FMin.HasReachedCallLimit()),
      EdmAboveMax(FMin.IsAboveMaxEdm()), HesseFailed(FMin.HesseFailed()),
      ErrorDef(FMin.Up()), NFcn(FMin.NFcn()), Edm(FMin.Edm()),
      GlobalCC(getGlobalCorrelations(FMin.UserState())) {}

std::ostream &operator<<(std::ostream &os, const MinuitResult &Result) {
  const FitResult &BaseResult = Result;
  os << BaseResult;

  os << "--------------MINUIT2 FIT INFOS----------------\n";
  if (!Result.IsValid)
    os << "    *** MINIMUM NOT VALID! ***\n";
  os << "Estimated distance to minimum: " << Result.Edm << "\n";
  if (Result.EdmAboveMax)
    os << "    *** EDM IS ABOVE MAXIMUM! ***\n";
  os << "Error definition: " << Result.ErrorDef << "\n";
  os << "Number of calls: " << Result.NFcn << "\n";
  if (Result.HasReachedCallLimit)
    os << "    *** LIMIT OF MAX CALLS REACHED! ***\n";
  if (!Result.HasValidParameters)
    os << "    *** NO VALID SET OF PARAMETERS! ***\n";
  if (!Result.HasValidCov)
    os << "    *** COVARIANCE MATRIX NOT VALID! ***\n";
  if (!Result.HasAccCov)
    os << "    *** COVARIANCE MATRIX NOT ACCURATE! ***\n";
  if (!Result.CovPosDef)
    os << "    *** COVARIANCE MATRIX NOT POSITIVE DEFINITE! ***\n";
  if (Result.HesseFailed)
    os << "    *** HESSE FAILED! ***\n";
  os << "-----------------------------------------------\n";
  return os;
}

void MinuitResult::write(std::string filename) const {
  std::ofstream ofs(filename);
  boost::archive::xml_oarchive oa(ofs);
  oa << boost::serialization::make_nvp("MinuitResult", *this);
}

MinuitResult load(std::string filename) {
  MinuitResult Result;
  std::ifstream ifs(filename);
  assert(ifs.good());
  boost::archive::xml_iarchive ia(ifs);
  ia >> boost::serialization::make_nvp("MinuitResult", Result);
  return Result;
}

} // namespace Minuit2
} // namespace Optimizer
} // namespace ComPWA
