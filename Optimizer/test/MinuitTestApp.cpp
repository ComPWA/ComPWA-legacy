// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE MinuitTests

#include <boost/test/unit_test.hpp>

#include "Core/Logging.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Estimator/Estimator.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

BOOST_AUTO_TEST_SUITE(OptimizerTests)

// this function just wraps a std::function and evaluates them with the
class FunctionRootEstimator : public ComPWA::Estimator::Estimator {
public:
  FunctionRootEstimator(std::function<double(double)> f,
                        const ComPWA::ParameterList &pars)
      : Function(f), Parameters(pars) {}

  double evaluate() const {
    return std::abs(Function(Parameters.doubleParameter(0)->value()));
  }

private:
  std::function<double(double)> Function;
  ComPWA::ParameterList Parameters;
};

void testFunctionRootFind(std::function<double(double)> f, double StartValue,
                          double ExpectedResult) {
  // linear function
  auto x = std::make_shared<ComPWA::FitParameter>("x", StartValue);
  x->fixParameter(false);
  ComPWA::ParameterList Params;
  Params.addParameter(x);

  auto LinearRootFind = std::make_shared<FunctionRootEstimator>(f, Params);

  auto Minimizer = std::make_shared<ComPWA::Optimizer::Minuit2::MinuitIF>(
      LinearRootFind, Params);
  auto Result = Minimizer->exec(Params);

  LOG(DEBUG) << Result->finalParameters().doubleParameter(0)->value();
  BOOST_CHECK(std::abs(Result->finalParameters().doubleParameter(0)->value() -
                       ExpectedResult) < 1e-4);
}

/// We just try to find the root of a function using minuit2
BOOST_AUTO_TEST_CASE(Minuit2OptimizationTest) {
  testFunctionRootFind([](double x) { return x - 2.0; }, 10.0, 2.0);
  testFunctionRootFind([](double x) { return x * x - 1.0; }, 4.0, 1.0);
  testFunctionRootFind([](double x) { return 3 * x * x - 2 * x - 4.0; }, -5.0,
                       -0.868517);
};
BOOST_AUTO_TEST_SUITE_END()
