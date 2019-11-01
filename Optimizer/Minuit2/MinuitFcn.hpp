// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef OPTIMIZER_MINUIT2_MINUITFCN_HPP_
#define OPTIMIZER_MINUIT2_MINUITFCN_HPP_

#include "Estimator/Estimator.hpp"

#include "Minuit2/FCNBase.h"

#include <iomanip>
#include <map>
#include <sstream>

namespace ROOT {
namespace Minuit2 {

///
/// \class MinuitFcn
/// Minuit2 function to be optimized based on the Minuit2 FcnBase. This class
/// uses the Estimator interface for the optimization.
///
class MinuitFcn : public FCNBase {

public:
  MinuitFcn(ComPWA::Estimator::Estimator<double> &estimator)
      : Estimator(estimator){};
  virtual ~MinuitFcn() = default;

  double operator()(const std::vector<double> &x) const {
    Estimator.updateParametersFrom(x);

    // Start timing
    std::chrono::steady_clock::time_point StartTime =
        std::chrono::steady_clock::now();
    double result = Estimator.evaluate();
    std::chrono::steady_clock::time_point EndTime =
        std::chrono::steady_clock::now();

    LOG(DEBUG) << "MinuitFcn: Estimator = " << std::setprecision(10) << result
               << std::setprecision(4) << " Time: "
               << std::chrono::duration_cast<std::chrono::milliseconds>(
                      EndTime - StartTime)
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
  };

  double Up() const {
    return 0.5; // TODO: Setter, LH 0.5, Chi2 1.
  };

private:
  ComPWA::Estimator::Estimator<double> &Estimator;
};

} // namespace Minuit2
} // namespace ROOT

#endif
