// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cassert>
#include <chrono>
#include <iomanip>

#include "Core/FitParameter.hpp"
#include "Core/Logging.hpp"
#include "Core/ParameterList.hpp"
#include "Estimator/Estimator.hpp"
#include "MinuitFcn.hpp"

using namespace ROOT::Minuit2;

MinuitFcn::MinuitFcn(std::shared_ptr<ComPWA::Estimator::Estimator> estimator,
                     ComPWA::ParameterList &parameters)
    : Estimator(estimator), Parameters(parameters) {
  if (0 == Estimator)
    throw std::runtime_error(
        "MinuitFcn::MinuitFcn() | Estimator is uninitialized!");
}

MinuitFcn::~MinuitFcn() {}

double MinuitFcn::operator()(const std::vector<double> &x) const {
  std::ostringstream paramOut;

  size_t pos = 0;
  for (auto p : Parameters.doubleParameters()) {
    if (p->isFixed()) {
      ++pos;
      continue;
    }
    paramOut << x[pos] << " "; // print only free parameters
    p->setValue(x[pos]);
    ++pos;
  }
  assert(x.size() == pos && "MinuitFcn::operator() | Number is (internal) "
                            "Minuit parameters and number of ComPWA "
                            "parameters does not match!");

  // Start timing
  std::chrono::steady_clock::time_point StartTime =
      std::chrono::steady_clock::now();
  double result = Estimator->evaluate();
  std::chrono::steady_clock::time_point EndTime =
      std::chrono::steady_clock::now();

  ;
  LOG(DEBUG) << "MinuitFcn: Estimator = " << std::setprecision(10) << result
             << std::setprecision(4) << " Time: "
             << std::chrono::duration_cast<std::chrono::milliseconds>(EndTime -
                                                                      StartTime)
                    .count()
             << "ms";
  LOG(DEBUG) << "Parameters: " << paramOut.str();

  return result;
}

double MinuitFcn::Up() const {
  return 0.5; // TODO: Setter, LH 0.5, Chi2 1.
}
