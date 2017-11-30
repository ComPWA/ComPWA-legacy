// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <vector>
#include <ctime>
#include <string>
#include <sstream>
#include <iostream>
#include <memory>

#include <boost/archive/xml_iarchive.hpp>

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MinosError.h"

#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"
#include "Core/FitResult.hpp"

using namespace ComPWA::Optimizer::Minuit2;

using namespace boost::log;
using namespace ROOT::Minuit2;

double shiftAngle(double v) {
  double originalVal = v;
  double val = originalVal;
  while (val > M_PI)
    val -= 2 * M_PI;
  while (val < (-1) * M_PI)
    val += 2 * M_PI;
  if (val != originalVal)
    LOG(info) << "shiftAngle() | Shifting parameter from " << originalVal
              << " to " << val << "!";
  return val;
}

MinuitIF::MinuitIF(std::shared_ptr<IEstimator> esti, ParameterList &par)
    : _myFcn(esti, par), estimator(esti), enableHesse(true), enableMinos(true) {

}

MinuitIF::~MinuitIF() {}

std::shared_ptr<ComPWA::FitResult> MinuitIF::exec(ParameterList &par) {
  LOG(debug) << "MinuitIF::exec() | Start";

  // Start timing
  clock_t begin = clock();

  par.RemoveDuplicates();

  LOG(debug) << "MinuitIF::exec() | Begin ParameterList::DeepCopy()";

  ParameterList initialParList;
  initialParList.DeepCopy(par);

  MnUserParameters upar;
  int freePars = 0;
  for (unsigned int i = 0; i < par.GetNDouble(); ++i) {
    std::shared_ptr<DoubleParameter> actPat = par.GetDoubleParameter(i);

    // If no error is set or error set to 0 we use a default error,
    // otherwise minuit treads this parameter as fixed
    if (!actPat->HasError() || actPat->GetError() <= 0)
      actPat->SetError(0.001);

    if (!actPat->IsFixed() &&
        actPat->name().find("phase") != actPat->name().npos)
      actPat->SetValue(shiftAngle(actPat->GetValue()));

    if (actPat->HasBounds()) {
      upar.Add(actPat->name(), actPat->GetValue(), actPat->GetError(),
               actPat->GetMinValue(), actPat->GetMaxValue());
    } else {
      upar.Add(actPat->name(), actPat->GetValue(), actPat->GetError());
    }

    if (!actPat->IsFixed())
      freePars++;
    if (actPat->IsFixed())
      upar.Fix(actPat->name());
  }

  LOG(info) << "MinuitIF::exec() | Number of parameters (free): "
            << par.GetNDouble() << " (" << freePars << ")";

  // use MnStrategy class to set all options for the fit
  MinuitStrategy strat(1); // using default strategy = 1 (medium)

  try { // try to read xml file for minuit setting
    std::ifstream ifs("MinuitStrategy.xml");
    boost::archive::xml_iarchive ia(ifs, boost::archive::no_header);
    ia >> BOOST_SERIALIZATION_NVP(strat);
    strat.init(); // update parameters of MnStrategy mother class (IMPORTANT!)
    ifs.close();
    LOG(debug) << "Minuit strategy parameters: from MinuitStrategy.xml";
  } catch (std::exception &ex) {
  }

  LOG(debug) << "Gradient number of steps: " << strat.GradientNCycles();
  LOG(debug) << "Gradient step tolerance: " << strat.GradientStepTolerance();
  LOG(debug) << "Gradient tolerance: " << strat.GradientTolerance();
  LOG(debug) << "Hesse number of steps: " << strat.HessianNCycles();
  LOG(debug) << "Hesse gradient number of steps: "
             << strat.HessianGradientNCycles();
  LOG(debug) << "Hesse step tolerance: " << strat.HessianStepTolerance();
  LOG(debug) << "Hesse G2 tolerance: " << strat.HessianG2Tolerance();

  // MIGRAD
  MnMigrad migrad(_myFcn, upar, strat);
  double maxfcn = 0.0;
  double tolerance = 0.1;

  /* From the MINUIT2 Documentation:
   * Minimize the function MnMigrad()(maxfcn, tolerance)
   *    @param maxfcn : max number of function calls (if = 0) default is
   *           used which is set to 200 + 100 * npar + 5 * npar**2
   *    @param tolerance : value used for terminating iteration procedure.
   *           For example, MIGRAD will stop iterating when edm (expected
   *           distance from minimum) will be: edm < tolerance * 10**-3
   *           Default value of tolerance used is 0.1
   */
  LOG(info) << "MinuitIF::exec() | Starting migrad: "
               "maxCalls="
            << maxfcn << " tolerance=" << tolerance;

  FunctionMinimum minMin = migrad(maxfcn, tolerance); //(maxfcn,tolerance)

  LOG(info) << "MinuitIF::exec() | Migrad finished! "
               "Minimum is valid = "
            << minMin.IsValid();

  // HESSE
  MnHesse hesse(strat);
  if (minMin.IsValid() && enableHesse) {
    LOG(info) << "MinuitIF::exec() | Starting hesse";
    hesse(_myFcn, minMin); // function minimum minMin is updated by hesse
    LOG(info) << "MinuitIF::exec() | Hesse finished";
  } else
    LOG(info) << "MinuitIF::exec() | Migrad failed to "
                 "find minimum! Skip hesse and minos!";

  LOG(info) << "MinuitIF::exec() | Minimization finished! "
               "LH = "
            << std::setprecision(10) << minMin.Fval();

  // MINOS
  MnMinos minos(_myFcn, minMin, strat);

  // save minimzed values
  MnUserParameterState minState = minMin.UserState();

  // ParameterList can be changed by minos. We have to do a deep copy here
  // to preserve the original parameters at the minimum.
  ParameterList finalParList;
  finalParList.DeepCopy(par);

  // We directly write the central values of the fit result to check if the
  // parameters change later on
  std::stringstream resultsOut;
  resultsOut << "Central values of floating paramters:" << std::endl;
  for (unsigned int i = 0; i < finalParList.GetNDouble(); ++i) {
    auto finalPar = finalParList.GetDoubleParameter(i);
    if (finalPar->IsFixed())
      continue;

    // central value
    double val = minState.Value(finalPar->name());

    // shift to [-pi;pi] if parameter is a phase
    if (finalPar->name().find("phase") != finalPar->name().npos)
      val = shiftAngle(val);
    finalPar->SetValue(val);

    resultsOut << finalPar->name() << " " << val << std::endl;
    if (finalPar->GetErrorType() == ErrorType::ASYM) {
      // Skip minos and fill symmetic errors
      if (!minMin.IsValid() || !enableMinos) {
        LOG(info) << "MinuitIF::exec() | Skip Minos "
                     "for parameter "
                  << i << "...";
        finalPar->SetError(minState.Error(finalPar->name()));
        continue;
      }
      // asymmetric errors -> run minos
      LOG(info) << "MinuitIF::exec() | Run minos "
                   "for parameter "
                << i << "...";
      MinosError err = minos.Minos(i);
      // lower = pair.first, upper= pair.second
      std::pair<double, double> assymErrors = err();
      finalPar->SetError(assymErrors.first, assymErrors.second);
    } else if (finalPar->GetErrorType() == ErrorType::SYM) {
      // symmetric errors -> migrad/hesse error
      finalPar->SetError(minState.Error(finalPar->name()));
    } else {
      throw std::runtime_error(
          "MinuitIF::exec() | Unknown error type of parameter: " +
          std::to_string((long long int)finalPar->GetErrorType()));
    }
  }

  // Update the original parameter list
  for (unsigned int i = 0; i < finalParList.GetNDouble(); ++i) {
    auto finalPar = finalParList.GetDoubleParameter(i);
    if (finalPar->IsFixed())
      continue;
    par.GetDoubleParameter(i)->UpdateParameter(finalPar);
  }

  LOG(debug) << "MinuitIF::exec() | " << resultsOut.str();

  double elapsed = double(clock() - begin) / CLOCKS_PER_SEC;

  // Create fit result
  std::shared_ptr<FitResult> result(new MinuitResult(estimator, minMin));
  result->SetFinalParameters(finalParList);
  result->SetInitialParameters(initialParList);
  result->SetTime(elapsed);

  // update parameters in amplitude
  //  Amplitude::UpdateAmpParameterList(estimator->getAmplitudes(),
  //  finalParList);

  return result;
}
