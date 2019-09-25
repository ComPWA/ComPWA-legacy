// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE GenevaTests

#include "Core/FitParameter.hpp"
#include "Core/Logging.hpp"
#include "Estimator/Estimator.hpp"
#include "Optimizer/Geneva/GenevaIF.hpp"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(OptimizerTests)

// this function just wraps a std::function and evaluates them with the
class FunctionRootEstimator1D : public ComPWA::Estimator::Estimator<double> {
public:
  FunctionRootEstimator1D(std::function<double(double)> f) : Function(f) {}

  double evaluate() noexcept { return std::abs(Function(x)); }

  void updateParametersFrom(const std::vector<double> &params) {
    x = params[0];
  }

  std::vector<ComPWA::Parameter> getParameters() const {
    return {ComPWA::Parameter{"x", x}};
  }

private:
  std::function<double(double)> Function;
  double x;
};

void testFunctionRootFind(std::function<double(double)> f, double StartValue,
                          std::vector<double> ExpectedResults) {
  // linear function
  ComPWA::FitParameter<double> x("x", StartValue, -20, 20, false);
  ComPWA::FitParameterList Params = {x};

  FunctionRootEstimator1D LinearRootFind(f);

  auto Minimizer = ComPWA::Optimizer::Geneva::GenevaIF();
  auto Result = Minimizer.optimize(LinearRootFind, Params);

  LOG(DEBUG) << Result;
  for (auto exp : ExpectedResults) {
    if (std::abs(Result.FinalParameters[0].Value - exp) < 1e-2) {
      BOOST_TEST(std::abs(Result.FinalParameters[0].Value - exp) < 1e-2);
      return;
    }
  }
  BOOST_TEST(false);
}

/// We just try to find the root of a function using geneva
BOOST_AUTO_TEST_CASE(GenevaOptimizationTest) {
  testFunctionRootFind([](double x) { return x - 2.0; }, 10.0, {2.0});
  testFunctionRootFind([](double x) { return x * x - 1.0; }, 4.0, {-1.0, 1.0});
  testFunctionRootFind([](double x) { return 3 * x * x - 2 * x - 4.0; }, -5.0,
                       {-0.868517});
};
BOOST_AUTO_TEST_SUITE_END()
