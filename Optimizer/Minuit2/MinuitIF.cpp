// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Core/FitParameter.hpp"
#include "Core/FitResult.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Utils.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnUserParameters.h"

#include "boost/archive/xml_iarchive.hpp"

using namespace ROOT::Minuit2;

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

MinuitIF::MinuitIF(std::shared_ptr<ComPWA::Estimator::Estimator> esti,
                   ParameterList &par, bool useHesse, bool useMinos)
    : Function(esti, par), Estimator(esti), FitStrategy(1), UseHesse(useHesse),
      UseMinos(useMinos) {
  // use MnStrategy class to set all options for the fit

  try { // try to read xml file for minuit setting
    std::ifstream ifs("MinuitStrategy.xml");
    boost::archive::xml_iarchive ia(ifs, boost::archive::no_header);
    ia >> BOOST_SERIALIZATION_NVP(FitStrategy);
    FitStrategy
        .init(); // update parameters of MnStrategy mother class (IMPORTANT!)
    ifs.close();
    LOG(DEBUG) << "Minuit strategy parameters: from MinuitStrategy.xml";

    LOG(DEBUG) << "Gradient number of steps: " << FitStrategy.GradientNCycles();
    LOG(DEBUG) << "Gradient step tolerance: "
               << FitStrategy.GradientStepTolerance();
    LOG(DEBUG) << "Gradient tolerance: " << FitStrategy.GradientTolerance();
    LOG(DEBUG) << "Hesse number of steps: " << FitStrategy.HessianNCycles();
    LOG(DEBUG) << "Hesse gradient number of steps: "
               << FitStrategy.HessianGradientNCycles();
    LOG(DEBUG) << "Hesse step tolerance: "
               << FitStrategy.HessianStepTolerance();
    LOG(DEBUG) << "Hesse G2 tolerance: " << FitStrategy.HessianG2Tolerance();
  } catch (std::exception &ex) {
  }
}

ComPWA::FitResult MinuitIF::execute(ParameterList &list) {
  LOG(DEBUG) << "MinuitIF::execute(): Start";

  std::chrono::steady_clock::time_point StartTime =
      std::chrono::steady_clock::now();

  LOG(DEBUG) << "MinuitIF::execute(): Begin ParameterList::DeepCopy()";

  ParameterList initialParList;
  initialParList.DeepCopy(list);

  MnUserParameters upar;
  int freePars = 0;
  for (auto actPat : list.doubleParameters()) {
    if (actPat->name() == "")
      throw BadParameter("MinuitIF::execute(): FitParameter without name in "
                         "list. Since FitParameter names are unique we stop "
                         "here.");
    // If no error is set or error set to 0 we use a default error,
    // otherwise minuit treads this parameter as fixed
    if (!actPat->hasError())
      actPat->setError(0.001);

    // Shift phase parameters to the range [-Pi,Pi]
    if (!actPat->isFixed() &&
        actPat->name().find("phase") != actPat->name().npos)
      actPat->setValue(
          ComPWA::Utils::shiftAngleToMinusPiPiWindow(actPat->value()));

    bool rt;
    if (actPat->hasBounds()) {
      rt = upar.Add(actPat->name(), actPat->value(), actPat->avgError(),
                    actPat->bounds().first, actPat->bounds().second);
      if (!rt)
        throw BadParameter(
            "MinuitIF::exec() | FitParameter " + actPat->name() +
            " can not be added to Minuit2 (internal) parameters.");
    } else {
      rt = upar.Add(actPat->name(), actPat->value(), actPat->avgError());
      if (!rt)
        throw BadParameter(
            "MinuitIF::exec() | FitParameter " + actPat->name() +
            " can not be added to Minuit2 (internal) parameters.");
    }
    if (!actPat->isFixed())
      freePars++;
    if (actPat->isFixed())
      upar.Fix(actPat->name());
  }

  LOG(INFO) << "MinuitIF::exec() | Number of parameters (free): "
            << list.numParameters() << " (" << freePars << ")";

  // MIGRAD
  MnMigrad migrad(Function, upar, FitStrategy);
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
  LOG(INFO) << "MinuitIF::exec() | Starting migrad: "
               "maxCalls="
            << maxfcn << " tolerance=" << tolerance;

  FunctionMinimum minMin = migrad(maxfcn, tolerance); //(maxfcn,tolerance)

  LOG(INFO) << "MinuitIF::exec() | Migrad finished! "
               "Minimum is valid = "
            << minMin.IsValid();

  // HESSE
  MnHesse hesse(FitStrategy);
  if (minMin.IsValid() && UseHesse) {
    LOG(INFO) << "MinuitIF::exec() | Starting hesse";
    hesse(Function, minMin); // function minimum minMin is updated by hesse
    LOG(INFO) << "MinuitIF::exec() | Hesse finished";
  } else
    LOG(INFO) << "MinuitIF::exec() | Migrad failed to "
                 "find minimum! Skip hesse and minos!";

  LOG(INFO) << "MinuitIF::exec() | Minimization finished! "
               "LH = "
            << std::setprecision(10) << minMin.Fval();

  std::chrono::steady_clock::time_point EndTime =
      std::chrono::steady_clock::now();

  // update the parameter list with the minuit2 parameters

  // calculate the errors using minos

  ParameterList FinalParameters = createFinalParameterList(list, minMin);
  FitResult result = createFitResult(FinalParameters, minMin);
  result.InitialParameters = initialParList;
  result.ElapsedTimeInSeconds =
      std::chrono::duration_cast<std::chrono::seconds>(EndTime - StartTime)
          .count();

  return result;
}

ParameterList MinuitIF::createFinalParameterList(
    ParameterList &list, const ROOT::Minuit2::FunctionMinimum &min) const {
  // ParameterList can be changed by minos. We have to do a deep copy here
  // to preserve the original parameters at the minimum.
  ParameterList finalParList;
  finalParList.DeepCopy(list);

  ROOT::Minuit2::MnUserParameterState minState = min.UserState();
  // MINOS
  MnMinos minos(Function, min, FitStrategy);

  // We directly write the central values of the fit result to check if the
  // parameters change later on
  std::stringstream resultsOut;
  resultsOut << "Central values of floating parameters:" << std::endl;
  size_t id = 0;
  for (auto finalPar : finalParList.doubleParameters()) {
    if (finalPar->isFixed())
      continue;
    // central value
    double val = minState.Value(finalPar->name());

    // shift to [-pi;pi] if parameter is a phase
    if (finalPar->name().find("phase") != finalPar->name().npos)
      val = ComPWA::Utils::shiftAngleToMinusPiPiWindow(val);
    finalPar->setValue(val);

    resultsOut << finalPar->name() << " " << val << std::endl;
    if (finalPar->errorType() == ErrorType::ASYM) {
      // Skip minos and fill symmetic errors
      if (!minState.IsValid() || !UseMinos) {
        LOG(INFO) << "MinuitIF::exec() | Skip Minos "
                     "for parameter "
                  << finalPar->name() << "...";
        finalPar->setError(minState.Error(finalPar->name()));
        continue;
      }
      // asymmetric errors -> run minos
      LOG(INFO) << "MinuitIF::exec() | Run minos "
                   "for parameter ["
                << id << "] " << finalPar->name() << "...";
      MinosError err = minos.Minos(id);
      // lower = pair.first, upper= pair.second
      std::pair<double, double> assymErrors = err();
      finalPar->setError(assymErrors.first, assymErrors.second);
    } else if (finalPar->errorType() == ErrorType::SYM) {
      // symmetric errors -> migrad/hesse error
      finalPar->setError(minState.Error(finalPar->name()));
    } else {
      throw std::runtime_error(
          "MinuitIF::exec() | Unknown error type of parameter: " +
          std::to_string((long long int)finalPar->errorType()));
    }
    id++;
  }

  LOG(DEBUG) << "MinuitIF::exec() | " << resultsOut.str();

  // Update the original parameter list
  for (size_t i = 0; i < finalParList.numParameters(); ++i) {
    auto finalPar = finalParList.doubleParameter(i);
    if (finalPar->isFixed())
      continue;
    list.doubleParameter(i)->updateParameter(finalPar);
  }
  return finalParList;
}

FitResult MinuitIF::createFitResult(const ParameterList &FinalParameterList,
                                    const ROOT::Minuit2::FunctionMinimum &min) {
  FitResult result;

  ROOT::Minuit2::MnUserParameterState minState = min.UserState();

  result.FinalParameters = FinalParameterList;

  size_t NumFreeParameter = minState.Parameters().Trafo().VariableParameters();

  if (minState.HasCovariance()) {
    ROOT::Minuit2::MnUserCovariance minuitCovMatrix = minState.Covariance();
    // Size of Minuit covariance vector is given by dim*(dim+1)/2.
    // dim is the dimension of the covariance matrix.
    // The dimension can therefore be calculated as
    // dim = -0.5+-0.5 sqrt(8*size+1)
    assert(minuitCovMatrix.Nrow() == NumFreeParameter);
    auto Cov = std::vector<std::vector<double>>(
        NumFreeParameter, std::vector<double>(NumFreeParameter));
    auto Corr = std::vector<std::vector<double>>(
        NumFreeParameter, std::vector<double>(NumFreeParameter));
    for (unsigned i = 0; i < NumFreeParameter; ++i)
      for (unsigned j = i; j < NumFreeParameter; ++j) {
        Cov.at(i).at(j) = minuitCovMatrix(j, i);
        Cov.at(j).at(i) = minuitCovMatrix(j, i); // fill lower half
      }
    for (unsigned i = 0; i < NumFreeParameter; ++i)
      for (unsigned j = i; j < NumFreeParameter; ++j) {
        Corr.at(i).at(j) =
            Cov.at(i).at(j) / sqrt(Cov.at(i).at(i) * Cov.at(j).at(j));
        Corr.at(j).at(i) = Corr.at(i).at(j); // fill lower half
      }

    result.MatrixProperties["Covariance"] = Cov;
    result.MatrixProperties["Correlation"] = Corr;
  } else {
    LOG(ERROR)
        << "MinuitIF::createFitResult: no valid correlation matrix available!";
  }
  if (minState.HasGlobalCC()) {
    result.DoubleListProperties["GlobalCC"] = minState.GlobalCC().GlobalCC();
  } else {
    LOG(ERROR)
        << "MinuitIF::createFitResult: no valid global correlation available!";
  }
  result.InitialEstimatorValue = -1.0;
  result.FinalEstimatorValue = minState.Fval();

  result.DoubleProperties["Edm"] = minState.Edm();
  result.IntProperties["IsValid"] = min.IsValid();
  result.IntProperties["HasPosDefCovariance"] = min.HasPosDefCovar();
  result.IntProperties["HasValidParameters"] = min.HasValidParameters();
  result.IntProperties["HasValidCovariance"] = min.HasValidCovariance();
  result.IntProperties["HasAccurateCovariance"] = min.HasAccurateCovar();
  result.IntProperties["HasReachedCallLimit"] = min.HasReachedCallLimit();
  result.IntProperties["EdmAboveMax"] = min.IsAboveMaxEdm();
  result.IntProperties["HesseFailed"] = min.HesseFailed();
  result.DoubleProperties["ErrorDef"] = min.Up();
  result.IntProperties["NFcn"] = min.NFcn();

  return result;
}

} // namespace Minuit2
} // namespace Optimizer
} // namespace ComPWA
