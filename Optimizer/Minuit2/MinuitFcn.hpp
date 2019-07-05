// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef OPTIMIZER_MINUIT2_MINUITFCN_HPP_
#define OPTIMIZER_MINUIT2_MINUITFCN_HPP_

#include "Estimator/Estimator.hpp"

#include "Minuit2/FCNBase.h"

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
  MinuitFcn(ComPWA::Estimator::Estimator<double> &estimator);
  virtual ~MinuitFcn() = default;

  double operator()(const std::vector<double> &x) const;

  double Up() const;

private:
  ComPWA::Estimator::Estimator<double> &Estimator;
};

} // namespace Minuit2
} // namespace ROOT

#endif
