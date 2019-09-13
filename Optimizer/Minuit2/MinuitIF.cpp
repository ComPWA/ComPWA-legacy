// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "MinuitIF.hpp"

#include "Core/Exceptions.hpp"
#include "Core/FitParameter.hpp"
#include "Core/Logging.hpp"
#include "Core/Utils.hpp"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>

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

void MinuitIF::useHesse(bool x) { UseHesse = x; }

void MinuitIF::useMinos(bool x) { UseMinos = x; }

MinuitIF::MinuitIF(bool UseHesse_, bool UseMinos_)
    : UseHesse(UseHesse_), UseMinos(UseMinos_) {}

MinuitResult MinuitIF::optimize(ComPWA::Estimator::Estimator<double> &Estimator,
                                ComPWA::FitParameterList InitialParameters) {
  LOG(DEBUG) << "MinuitIF::optimize() | Start";
  std::chrono::steady_clock::time_point StartTime =
      std::chrono::steady_clock::now();

  double InitialEstimatorValue(Estimator.evaluate());

  MnUserParameters upar;
  size_t FreeParameters(0);
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

  // use MnStrategy class to set all options for the fit
  MnStrategy strat(1); // using default strategy = 1 (medium)

  try { // try to read xml file for minuit setting
    std::ifstream ifs("MinuitStrategy.xml");
    boost::archive::xml_iarchive ia(ifs, boost::archive::no_header);
    ia >> BOOST_SERIALIZATION_NVP(strat);
    ifs.close();
    LOG(DEBUG) << "Minuit strategy parameters: from MinuitStrategy.xml";
  } catch (std::exception &ex) {
  }

  LOG(DEBUG) << "Gradient number of steps: " << strat.GradientNCycles();
  LOG(DEBUG) << "Gradient step tolerance: " << strat.GradientStepTolerance();
  LOG(DEBUG) << "Gradient tolerance: " << strat.GradientTolerance();
  LOG(DEBUG) << "Hesse number of steps: " << strat.HessianNCycles();
  LOG(DEBUG) << "Hesse gradient number of steps: "
             << strat.HessianGradientNCycles();
  LOG(DEBUG) << "Hesse step tolerance: " << strat.HessianStepTolerance();
  LOG(DEBUG) << "Hesse G2 tolerance: " << strat.HessianG2Tolerance();

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

FitParameterList MinuitIF::getFinalParameters(
    const ROOT::Minuit2::MnUserParameterState &minState,
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

std::vector<std::vector<double>> MinuitIF::getCovarianceMatrix(
    const ROOT::Minuit2::MnUserParameterState &minState) {
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

} // namespace Minuit2
} // namespace Optimizer
} // namespace ComPWA

template <typename Archive>
void boost::serialization::serialize(Archive &ar, ROOT::Minuit2::MnStrategy &s,
                                     const unsigned int version) {
  //    ar & BOOST_SERIALIZATION_NVP(fStrategy);
  unsigned int fGradNCyc(s.GradientNCycles());
  double fGradTlrStp(s.GradientStepTolerance());
  double fGradTlr(s.GradientTolerance());
  unsigned int fHessNCyc(s.HessianNCycles());
  double fHessTlrStp(s.HessianStepTolerance());
  double fHessTlrG2(s.HessianG2Tolerance());
  unsigned int fHessGradNCyc(s.HessianGradientNCycles());

  ar &BOOST_SERIALIZATION_NVP(fGradNCyc);
  ar &BOOST_SERIALIZATION_NVP(fGradTlrStp);
  ar &BOOST_SERIALIZATION_NVP(fGradTlr);
  ar &BOOST_SERIALIZATION_NVP(fHessNCyc);
  ar &BOOST_SERIALIZATION_NVP(fHessTlrStp);
  ar &BOOST_SERIALIZATION_NVP(fHessTlrG2);
  ar &BOOST_SERIALIZATION_NVP(fHessGradNCyc);

  s.SetGradientNCycles(fGradNCyc);
  s.SetGradientStepTolerance(fGradTlrStp);
  s.SetGradientTolerance(fGradTlr);
  s.SetHessianNCycles(fHessNCyc);
  s.SetHessianStepTolerance(fHessTlrStp);
  s.SetHessianG2Tolerance(fHessTlrG2);
  s.SetHessianGradientNCycles(fHessGradNCyc);
}
