// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE MinuitTests

#include <boost/test/unit_test.hpp>

#include "Core/Logging.hpp"
#include "Estimator/Estimator.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

BOOST_AUTO_TEST_SUITE(OptimizerTests)

// this function just wraps a std::function and evaluates them with the
class FunctionRootEstimator : public ComPWA::Estimator::Estimator<double> {
public:
  FunctionRootEstimator(std::function<double(double)> f,
                        const std::vector<double> &pars)
      : Function(f), Parameters(pars) {}

  double evaluate() noexcept { return std::abs(Function(Parameters[0])); }

  void updateParametersFrom(const std::vector<double> &params) {
    Parameters = params;
  }
  std::vector<ComPWA::Parameter> getParameters() const {
    std::vector<ComPWA::Parameter> Params;
    size_t counter(0);
    for (auto x : Parameters) {
      Params.push_back(ComPWA::Parameter{"p" + counter, x});
    }
    return Params;
  }

private:
  std::function<double(double)> Function;
  std::vector<double> Parameters;
};

void testFunctionRootFind(std::function<double(double)> f, double StartValue,
                          double ExpectedResult) {
  // linear function
  auto x = ComPWA::FitParameter<double>("x", StartValue, false);

  auto LinearRootFind = FunctionRootEstimator(f, {StartValue});
  ComPWA::FitParameterList Params;
  Params.push_back(x);

  auto Minimizer = ComPWA::Optimizer::Minuit2::MinuitIF();
  auto Result = Minimizer.optimize(LinearRootFind, Params);

  LOG(DEBUG) << "root value: " << Result.FinalParameters[0].Value;
  BOOST_CHECK(std::abs(Result.FinalParameters[0].Value - ExpectedResult) <
              1e-4);
}

/// We just try to find the root of a function using minuit2
BOOST_AUTO_TEST_CASE(Minuit2OptimizationTest) {
  testFunctionRootFind([](double x) { return x - 2.0; }, 10.0, 2.0);
  testFunctionRootFind([](double x) { return x * x - 1.0; }, 4.0, 1.0);
  testFunctionRootFind([](double x) { return 3 * x * x - 2 * x - 4.0; }, -5.0,
                       -0.868517);
};
BOOST_AUTO_TEST_SUITE_END()
