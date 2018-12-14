// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <ctime>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnUserParameters.h"

#include "Core/FitParameter.hpp"
#include "Core/FitResult.hpp"
#include "Core/ParameterList.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

using namespace ComPWA::Optimizer::Minuit2;

using namespace ROOT::Minuit2;

double shiftAngle(double v) {
  double originalVal = v;
  double val = originalVal;
  while (val > M_PI)
    val -= 2 * M_PI;
  while (val < (-1) * M_PI)
    val += 2 * M_PI;
  if (val != originalVal)
    LOG(INFO) << "shiftAngle() | Shifting parameter from " << originalVal
              << " to " << val << "!";
  return val;
}

MinuitIF::MinuitIF(std::shared_ptr<ComPWA::Estimator::Estimator> esti,
                   ParameterList &par)
    : Function(esti, par), Estimator(esti), UseHesse(true), UseMinos(true) {}

MinuitIF::~MinuitIF() {}

std::shared_ptr<ComPWA::FitResult> MinuitIF::exec(ParameterList &list) {
  LOG(DEBUG) << "MinuitIF::exec() | Start";

  // Start timing
  clock_t begin = clock();

  LOG(DEBUG) << "MinuitIF::exec() | Begin ParameterList::DeepCopy()";

  ParameterList initialParList;
  initialParList.DeepCopy(list);

  MnUserParameters upar;
  int freePars = 0;
  for (auto actPat : list.doubleParameters()) {
    if (actPat->name() == "")
      throw BadParameter("MinuitIF::exec() | FitParameter without name in "
                         "list. Since FitParameter names are unique we stop "
                         "here.");
    // If no error is set or error set to 0 we use a default error,
    // otherwise minuit treads this parameter as fixed
    if (!actPat->hasError())
      actPat->setError(0.001);

    // Shift phase parameters to the range [-Pi,Pi]
    if (!actPat->isFixed() &&
        actPat->name().find("phase") != actPat->name().npos)
      actPat->setValue(shiftAngle(actPat->value()));

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

  // use MnStrategy class to set all options for the fit
  MinuitStrategy strat(1); // using default strategy = 1 (medium)

  try { // try to read xml file for minuit setting
    std::ifstream ifs("MinuitStrategy.xml");
    boost::archive::xml_iarchive ia(ifs, boost::archive::no_header);
    ia >> BOOST_SERIALIZATION_NVP(strat);
    strat.init(); // update parameters of MnStrategy mother class (IMPORTANT!)
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
  LOG(INFO) << "MinuitIF::exec() | Starting migrad: "
               "maxCalls="
            << maxfcn << " tolerance=" << tolerance;

  FunctionMinimum minMin = migrad(maxfcn, tolerance); //(maxfcn,tolerance)

  LOG(INFO) << "MinuitIF::exec() | Migrad finished! "
               "Minimum is valid = "
            << minMin.IsValid();

  // HESSE
  MnHesse hesse(strat);
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

  // MINOS
  MnMinos minos(Function, minMin, strat);

  // save minimzed values
  MnUserParameterState minState = minMin.UserState();

  // ParameterList can be changed by minos. We have to do a deep copy here
  // to preserve the original parameters at the minimum.
  ParameterList finalParList;
  finalParList.DeepCopy(list);

  // We directly write the central values of the fit result to check if the
  // parameters change later on
  std::stringstream resultsOut;
  resultsOut << "Central values of floating paramters:" << std::endl;
  size_t id = 0;
  for (auto finalPar : finalParList.doubleParameters()) {
    if (finalPar->isFixed())
      continue;
    // central value
    double val = minState.Value(finalPar->name());

    // shift to [-pi;pi] if parameter is a phase
    if (finalPar->name().find("phase") != finalPar->name().npos)
      val = shiftAngle(val);
    finalPar->setValue(val);

    resultsOut << finalPar->name() << " " << val << std::endl;
    if (finalPar->errorType() == ErrorType::ASYM) {
      // Skip minos and fill symmetic errors
      if (!minMin.IsValid() || !UseMinos) {
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

  // Update the original parameter list
  for (size_t i = 0; i < finalParList.numParameters(); ++i) {
    auto finalPar = finalParList.doubleParameter(i);
    if (finalPar->isFixed())
      continue;
    list.doubleParameter(i)->updateParameter(finalPar);
  }

  LOG(DEBUG) << "MinuitIF::exec() | " << resultsOut.str();

  double elapsed = double(clock() - begin) / CLOCKS_PER_SEC;

  // Create fit result
  std::shared_ptr<FitResult> result(new MinuitResult(Estimator, minMin));
  result->setFinalParameters(finalParList);
  result->setInitialParameters(initialParList);
  result->setTime(elapsed);

  // update parameters in amplitude
  //  Amplitude::UpdateAmpParameterList(estimator->getAmplitudes(),
  //  finalParList);

  return result;
}
