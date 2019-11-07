// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "MinuitIF.hpp"
#include "MinuitFcn.hpp"

#include "Core/Exceptions.hpp"
#include "Core/FitParameter.hpp"
#include "Core/Logging.hpp"
#include "Core/Utils.hpp"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnUserParameters.h"

#include <chrono>
#include <iomanip>
#include <string>
#include <vector>

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

using namespace ROOT::Minuit2;

FitParameterList
getFinalParameters(const ROOT::Minuit2::MnUserParameterState &minState,
                   FitParameterList InitialParameters) {
  FitParameterList FinalParameters(InitialParameters);

  for (auto &FinalPar : FinalParameters) {
    if (FinalPar.IsFixed)
      continue;

    // central value
    double val = minState.Value(FinalPar.Name);
    double err = minState.Error(FinalPar.Name);

    // shift to [-pi;pi] if parameter is a phase
    if (FinalPar.Name.find("phase") != FinalPar.Name.npos)
      val = ComPWA::Utils::shiftAngle(val);
    FinalPar.Value = val;
    FinalPar.Error = std::make_pair(err, err);
  }
  return FinalParameters;
}

std::vector<std::vector<double>>
getCovarianceMatrix(const ROOT::Minuit2::MnUserParameterState &minState) {
  std::vector<std::vector<double>> CovarianceMatrix;

  if (minState.HasCovariance()) {
    auto NumFreeParameter = minState.Parameters().Trafo().VariableParameters();
    ROOT::Minuit2::MnUserCovariance minuitCovMatrix = minState.Covariance();
    // Size of Minuit covariance vector is given by dim*(dim+1)/2.
    // dim is the dimension of the covariance matrix.
    // The dimension can therefore be calculated as
    // dim = -0.5+-0.5 sqrt(8*size+1)
    assert(minuitCovMatrix.Nrow() == NumFreeParameter);
    CovarianceMatrix = std::vector<std::vector<double>>(
        NumFreeParameter, std::vector<double>(NumFreeParameter));
    for (unsigned i = 0; i < NumFreeParameter; ++i)
      for (unsigned j = i; j < NumFreeParameter; ++j) {
        CovarianceMatrix.at(i).at(j) = minuitCovMatrix(j, i);
        CovarianceMatrix.at(j).at(i) = minuitCovMatrix(j, i); // fill lower half
      }

  } else {
    LOG(ERROR)
        << "MinuitIF::createResult(): no valid covariance matrix available!";
  }
  return CovarianceMatrix;
}

void MinuitIF::setStrategy(std::string strategy) {
  MnStrategy strat; // using default strategy = 1 (medium)
  if (strategy == "low")
    strat.SetLowStrategy();
  else if (strategy == "medium")
    strat.SetMediumStrategy();
  else if (strategy == "high")
    strat.SetHighStrategy();
  else
    LOG(INFO) << "MinuitIF::setStrategy() | Minuit strategy must be "
                 "set to 'low', 'medium' or 'high'";

  GradientNCycles = strat.GradientNCycles();
  GradientStepTolerance = strat.GradientStepTolerance();
  GradientTolerance = strat.GradientTolerance();
  HessianNCycles = strat.HessianNCycles();
  HessianGradientNCycles = strat.HessianGradientNCycles();
  HessianStepTolerance = strat.HessianStepTolerance();
  HessianG2Tolerance = strat.HessianG2Tolerance();
}

std::string MinuitIF::checkStrategy() {
  MnStrategy strat;

  strat.SetLowStrategy();
  if (strat.GradientNCycles() == GradientNCycles &&
      strat.GradientStepTolerance() == GradientStepTolerance &&
      strat.GradientTolerance() == GradientTolerance &&
      strat.HessianNCycles() == HessianNCycles &&
      strat.HessianGradientNCycles() == HessianGradientNCycles &&
      strat.HessianStepTolerance() == HessianStepTolerance &&
      strat.HessianG2Tolerance() == HessianG2Tolerance)
    return "low";

  strat.SetMediumStrategy();
  if (strat.GradientNCycles() == GradientNCycles &&
      strat.GradientStepTolerance() == GradientStepTolerance &&
      strat.GradientTolerance() == GradientTolerance &&
      strat.HessianNCycles() == HessianNCycles &&
      strat.HessianGradientNCycles() == HessianGradientNCycles &&
      strat.HessianStepTolerance() == HessianStepTolerance &&
      strat.HessianG2Tolerance() == HessianG2Tolerance)
    return "medium";

  strat.SetHighStrategy();
  if (strat.GradientNCycles() == GradientNCycles &&
      strat.GradientStepTolerance() == GradientStepTolerance &&
      strat.GradientTolerance() == GradientTolerance &&
      strat.HessianNCycles() == HessianNCycles &&
      strat.HessianGradientNCycles() == HessianGradientNCycles &&
      strat.HessianStepTolerance() == HessianStepTolerance &&
      strat.HessianG2Tolerance() == HessianG2Tolerance)
    return "high";

  return "custom";
}

MinuitResult MinuitIF::optimize(ComPWA::Estimator::Estimator<double> &Estimator,
                                ComPWA::FitParameterList InitialParameters) {
  LOG(DEBUG) << "MinuitIF::optimize() | Start";
  std::chrono::steady_clock::time_point StartTime =
      std::chrono::steady_clock::now();

  double InitialEstimatorValue(Estimator.evaluate());

  MnUserParameters upar;
  unsigned int FreeParameters = 0;
  for (auto Param : InitialParameters) {
    if (Param.Name == "")
      throw BadParameter("MinuitIF::optimize() | FitParameter without name in "
                         "list. Since FitParameter names are unique we stop "
                         "here.");
    // If no error is set or error set to 0 we use a default error,
    // otherwise minuit treats this parameter as fixed
    double Error(0.5 * (Param.Error.first + Param.Error.second));
    if (Error == 0.0) {
      if (Param.Value != 0.0)
        Error = Param.Value * 0.1;
      else
        Error = 0.01;
    }

    if (Param.HasBounds) {
      bool rt = upar.Add(Param.Name, Param.Value, Error, Param.Bounds.first,
                         Param.Bounds.second);
      if (!rt)
        throw BadParameter(
            "MinuitIF::optimize() | FitParameter " + Param.Name +
            " can not be added to Minuit2 (internal) parameters.");
    } else {
      bool rt = upar.Add(Param.Name, Param.Value, Error);
      if (!rt)
        throw BadParameter(
            "MinuitIF::optimize() | FitParameter " + Param.Name +
            " can not be added to Minuit2 (internal) parameters.");
    }
    if (!Param.IsFixed)
      FreeParameters++;
    if (Param.IsFixed)
      upar.Fix(Param.Name);
  }

  ROOT::Minuit2::MinuitFcn Function(Estimator);

  LOG(INFO) << "MinuitIF::optimize() | Number of parameters (free): "
            << InitialParameters.size() << " (" << FreeParameters << ")";

  // Configure minimization strategy of Minuit2
  MnStrategy strat; // default strategy = 1 (medium)

  strat.SetGradientNCycles(GradientNCycles);
  strat.SetGradientStepTolerance(GradientStepTolerance);
  strat.SetGradientTolerance(GradientTolerance);
  strat.SetHessianNCycles(HessianNCycles);
  strat.SetHessianGradientNCycles(HessianGradientNCycles);
  strat.SetHessianStepTolerance(HessianStepTolerance);
  strat.SetHessianG2Tolerance(HessianG2Tolerance);

  LOG(INFO) << "Minuit2 strategy: " << checkStrategy()
            << "\n       GradientNCycles: " << GradientNCycles
            << "\n       GradientStepTolerance: " << GradientStepTolerance
            << "\n       GradientTolerance: " << GradientTolerance
            << "\n       HessianNCycles: " << HessianNCycles
            << "\n       HessianGradientNCycles: " << HessianGradientNCycles
            << "\n       HessianStepTolerance: " << HessianStepTolerance
            << "\n       HessianG2Tolerance: " << HessianG2Tolerance;

  // MIGRAD
  MnMigrad migrad(Function, upar, strat);
  double maxfcn = 0.0;
  double tolerance = 0.1;

  // From the MINUIT2 Documentation:
  // Minimize the function MnMigrad()(maxfcn, tolerance)
  // \param maxfcn : max number of function calls (if = 0) default is
  //           used which is set to 200 + 100 * npar + 5 * npar**2
  // \param tolerance : value used for terminating iteration procedure.
  //           For example, MIGRAD will stop iterating when edm (expected
  //           distance from minimum) will be: edm < tolerance * 10**-3
  //           Default value of tolerance used is 0.1
  LOG(INFO) << "MinuitIF::optimize() | Starting migrad: "
               "maxCalls="
            << maxfcn << " tolerance=" << tolerance;

  FunctionMinimum minMin = migrad(maxfcn, tolerance); //(maxfcn,tolerance)

  LOG(INFO) << "MinuitIF::optimize() | Migrad finished! "
               "Minimum is valid = "
            << minMin.IsValid();

  // HESSE
  MnHesse hesse(strat);
  if (minMin.IsValid() && UseHesse) {
    LOG(INFO) << "MinuitIF::optimize() | Starting hesse";
    hesse(Function, minMin); // function minimum minMin is updated by hesse
    LOG(INFO) << "MinuitIF::optimize() | Hesse finished";
  } else
    LOG(INFO) << "MinuitIF::optimize() | Migrad failed to "
                 "find minimum! Skip hesse and minos!";

  LOG(INFO) << "MinuitIF::optimize() | Minimization finished! "
               "LH = "
            << std::setprecision(10) << minMin.Fval();

  std::chrono::steady_clock::time_point EndTime =
      std::chrono::steady_clock::now();

  // Create fit result, before changes to parameters in
  // FunctionMinimum might occur
  auto MinState = minMin.UserState();
  FitResult BaseResult{
      InitialParameters,
      getFinalParameters(MinState, InitialParameters),
      FreeParameters,
      minMin.IsValid(),
      InitialEstimatorValue,
      MinState.Fval(),
      std::chrono::duration_cast<std::chrono::seconds>(EndTime - StartTime),
      getCovarianceMatrix(MinState)};
  MinuitResult Result(BaseResult, minMin);

  // In case Minos should be used, recalculate the errors
  if (minMin.IsValid() && UseMinos) {
    // MINOS
    MnMinos minos(Function, minMin, strat);
    size_t id = 0;
    for (auto &FinalPar : Result.FinalParameters) {
      if (FinalPar.IsFixed)
        continue;

      // asymmetric errors -> run minos
      LOG(INFO) << "MinuitIF::optimize() | Run minos "
                   "for parameter ["
                << id << "] " << FinalPar.Name << "...";
      MinosError err = minos.Minos(id);
      FinalPar.Error = err();
      ++id;
    }
  }

  return Result;
}

} // namespace Minuit2
} // namespace Optimizer
} // namespace ComPWA
