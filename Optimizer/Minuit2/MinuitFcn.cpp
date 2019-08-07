// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <chrono>
#include <iomanip>

#include "Core/FitParameter.hpp"
#include "Core/Logging.hpp"
#include "MinuitFcn.hpp"

using namespace ROOT::Minuit2;

MinuitFcn::MinuitFcn(ComPWA::Estimator::Estimator<double> &estimator)
    : Estimator(estimator) {}

double MinuitFcn::operator()(const std::vector<double> &x) const {
  Estimator.updateParametersFrom(x);

  // Start timing
  std::chrono::steady_clock::time_point StartTime =
      std::chrono::steady_clock::now();
  double result = Estimator.evaluate();
  std::chrono::steady_clock::time_point EndTime =
      std::chrono::steady_clock::now();

  LOG(DEBUG) << "MinuitFcn: Estimator = " << std::setprecision(10) << result
             << std::setprecision(4) << " Time: "
             << std::chrono::duration_cast<std::chrono::milliseconds>(EndTime -
                                                                      StartTime)
                    .count()
             << "ms";
  LOG(DEBUG) << "Parameters: " << [&]() {
    std::ostringstream params;
    for (auto var : x) {
      params << var << " ";
    }
    return params.str();
  }();

  return result;
}

double MinuitFcn::Up() const {
  return 0.5; // TODO: Setter, LH 0.5, Chi2 1.
}
