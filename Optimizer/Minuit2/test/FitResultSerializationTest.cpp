// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE MinuitTests

#include <boost/test/unit_test.hpp>

#include <random>

#include "Optimizer/Minuit2/MinuitResult.hpp"

BOOST_AUTO_TEST_SUITE(OptimizerTests)

std::pair<ComPWA::FitParameterList, ComPWA::FitParameterList>
createRandomFitParameters() {
  std::random_device Random;
  std::mt19937 RandomGenerator(Random());
  std::uniform_int_distribution<size_t> SizeDistribution(1, 100);
  std::uniform_int_distribution<size_t> BoolDistribution(0, 1);
  std::uniform_real_distribution<> ValueDistribution(-100.0, 100.0);
  size_t NumVars = SizeDistribution(RandomGenerator);
  ComPWA::FitParameterList InitialParameters;
  ComPWA::FitParameterList FinalParameters;
  for (size_t i = 0; i < NumVars; ++i) {
    std::stringstream name;
    name << "Var" << i + 1;
    auto Var = ComPWA::FitParameter<double>(name.str(),
                                            ValueDistribution(RandomGenerator));
    Var.HasBounds = BoolDistribution(RandomGenerator);
    if (Var.HasBounds) {
      double range = std::abs(ValueDistribution(RandomGenerator));
      Var.Bounds = std::make_pair(Var.Value - range, Var.Value + range);
    }
    Var.IsFixed = BoolDistribution(RandomGenerator);
    InitialParameters.push_back(Var);
  }
  for (ComPWA::FitParameter<double> x : InitialParameters) {
    if (!x.IsFixed) {
      if (x.HasBounds) {
        std::uniform_real_distribution<> BoundedValueDistribution(
            x.Bounds.first, x.Bounds.second);
        x.Value = BoundedValueDistribution(RandomGenerator);
      } else {
        x.Value = ValueDistribution(RandomGenerator);
      }
    }
    FinalParameters.push_back(x);
  }
  return std::make_pair(InitialParameters, FinalParameters);
}

ComPWA::Optimizer::Minuit2::MinuitResult createRandomMinuitResult() {
  std::random_device Random;
  std::mt19937 RandomGenerator(Random());
  std::uniform_int_distribution<size_t> SizeDistribution(1, 100);
  std::uniform_real_distribution<> ValueDistribution(-100.0, 100.0);
  ComPWA::Optimizer::Minuit2::MinuitResult Result;
  auto Params = createRandomFitParameters();
  Result.InitialEstimatorValue = ValueDistribution(RandomGenerator);
  Result.FinalEstimatorValue = ValueDistribution(RandomGenerator);
  Result.InitialParameters = Params.first;
  Result.FinalParameters = Params.second;
  Result.NFcn = SizeDistribution(RandomGenerator);
  Result.FitDuration = std::chrono::seconds(SizeDistribution(RandomGenerator));
  Result.CovarianceMatrix = {
      {ValueDistribution(RandomGenerator), ValueDistribution(RandomGenerator),
       ValueDistribution(RandomGenerator)},
      {ValueDistribution(RandomGenerator), ValueDistribution(RandomGenerator),
       ValueDistribution(RandomGenerator)},
      {ValueDistribution(RandomGenerator), ValueDistribution(RandomGenerator),
       ValueDistribution(RandomGenerator)}};

  return Result;
}

BOOST_AUTO_TEST_CASE(FitResultSerializationTest) {
  for (size_t i = 0; i < 100; ++i) {
    auto Result = createRandomMinuitResult();
    Result.write("test.xml");

    auto NewResult = ComPWA::Optimizer::Minuit2::load("test.xml");

    BOOST_CHECK(std::abs(Result.InitialEstimatorValue -
                         NewResult.InitialEstimatorValue) < 1e-6);
    BOOST_CHECK(std::abs(Result.FinalEstimatorValue -
                         NewResult.FinalEstimatorValue) < 1e-6);
    BOOST_CHECK(Result.InitialParameters.size() ==
                NewResult.InitialParameters.size());
    for (size_t i = 0; i < Result.InitialParameters.size(); ++i) {
      if (Result.InitialParameters.at(i).HasBounds) {
        BOOST_CHECK(std::abs(Result.InitialParameters.at(i).Bounds.first -
                             NewResult.InitialParameters.at(i).Bounds.first) <
                    1e-6);
        BOOST_CHECK(std::abs(Result.InitialParameters.at(i).Bounds.second -
                             NewResult.InitialParameters.at(i).Bounds.second) <
                    1e-6);
      }
      BOOST_CHECK(Result.InitialParameters.at(i).Name ==
                  NewResult.InitialParameters.at(i).Name);
      BOOST_CHECK(std::abs(Result.InitialParameters.at(i).Value -
                           NewResult.InitialParameters.at(i).Value) < 1e-6);
    }
    BOOST_CHECK(Result.FinalParameters.size() ==
                NewResult.FinalParameters.size());
    for (size_t i = 0; i < Result.FinalParameters.size(); ++i) {
      BOOST_CHECK(Result.FinalParameters.at(i).Name ==
                  NewResult.FinalParameters.at(i).Name);
      BOOST_CHECK(std::abs(Result.FinalParameters.at(i).Value -
                           NewResult.FinalParameters.at(i).Value) < 1e-6);
    }
    BOOST_CHECK(Result.NFcn == NewResult.NFcn);
    BOOST_CHECK(Result.FitDuration == NewResult.FitDuration);

    BOOST_CHECK(Result.CovarianceMatrix.size() ==
                NewResult.CovarianceMatrix.size());
    for (size_t i = 0; i < Result.CovarianceMatrix.size(); ++i) {
      for (size_t j = 0; j < Result.CovarianceMatrix.at(i).size(); ++j) {
        BOOST_CHECK(std::abs(Result.CovarianceMatrix.at(i).at(j) -
                             NewResult.CovarianceMatrix.at(i).at(j)) < 1e-6);
      }
    }
  }
};
BOOST_AUTO_TEST_SUITE_END()
