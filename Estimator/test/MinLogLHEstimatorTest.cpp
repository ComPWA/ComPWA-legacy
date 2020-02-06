// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Estimator_MinLogLHEstimatorTest

#include <chrono>
#include <cmath>
#include <random>

#include <boost/test/unit_test.hpp>

#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/FunctionTree/TreeNode.hpp"
#include "Data/DataSet.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"
#include "Tools/Integration.hpp"

// reference gaussian intensity
class Gaussian : public ComPWA::Intensity {
public:
  Gaussian(double mean, double width)
      : Mean(ComPWA::Parameter{"mean", mean}),
        Width(ComPWA::Parameter{"width", width}),
        Strength(ComPWA::Parameter{"strength", 1.0}) {}

  std::vector<double> evaluate(const ComPWA::DataMap &data) noexcept {
    auto const &xvals = data.at("x");
    std::vector<double> result(xvals.size());
    std::transform(xvals.begin(), xvals.end(), result.begin(), [&](double x) {
      return Strength.Value * std::exp(-0.5 * std::pow(x - Mean.Value, 2) /
                                       std::pow(Width.Value, 2));
    });
    return result;
  }

  void updateParametersFrom(const std::vector<double> &params) {
    if (params.size() != 3)
      throw std::runtime_error(
          "MinLogLHEstimatorTest_Gaussian::updateParametersFrom(): Parameter "
          "list size is incorrect!");
    Mean.Value = params[0];
    Width.Value = params[1];
    Strength.Value = params[2];
  }
  std::vector<ComPWA::Parameter> getParameters() const {
    return {Mean, Width, Strength};
  }

private:
  ComPWA::Parameter Mean;
  ComPWA::Parameter Width;
  ComPWA::Parameter Strength;
};

std::shared_ptr<ComPWA::FunctionTree::TreeNode>
createFunctionTree(std::shared_ptr<ComPWA::FunctionTree::FitParameter> mean,
                   std::shared_ptr<ComPWA::FunctionTree::FitParameter> width,
                   std::shared_ptr<ComPWA::FunctionTree::FitParameter> strength,
                   const ComPWA::FunctionTree::ParameterList &DataSample) {

  size_t n = DataSample.mDoubleValue(0)->values().size();

  using namespace ComPWA::FunctionTree;
  auto tr = std::make_shared<TreeNode>(
      MDouble("", n), std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  tr->addNode(createLeaf(strength));
  auto Exp = std::make_shared<TreeNode>(
      MDouble("", n),
      std::make_shared<ComPWA::FunctionTree::Exp>(ParType::MDOUBLE));
  tr->addNode(Exp);
  auto Exponent = std::make_shared<TreeNode>(
      MDouble("", n), std::make_shared<MultAll>(ParType::MDOUBLE));
  Exp->addNode(Exponent);
  Exponent->addNode(createLeaf(-0.5));

  auto Nominator = std::make_shared<TreeNode>(
      MDouble("", n), std::make_shared<Pow>(ParType::MDOUBLE, 2));
  Exponent->addNode(Nominator);
  auto Diff =
      std::make_shared<TreeNode>(std::make_shared<AddAll>(ParType::MDOUBLE));
  Nominator->addNode(Diff);
  Diff->addNode(createLeaf(DataSample.mDoubleValue(0)));
  auto Negate = std::make_shared<TreeNode>(
      std::shared_ptr<Parameter>(new Value<double>("negMean")),
      std::make_shared<MultAll>(ParType::DOUBLE));
  Diff->addNode(Negate);
  Negate->addNodes({createLeaf(-1.0), createLeaf(mean)});

  auto Inverse = std::make_shared<TreeNode>(
      std::shared_ptr<Parameter>(new Value<double>("invdenom")),
      std::make_shared<ComPWA::FunctionTree::Inverse>(ParType::DOUBLE));
  Exponent->addNode(Inverse);

  auto Denominator = std::make_shared<TreeNode>(
      std::shared_ptr<Parameter>(new Value<double>("denom")),
      std::make_shared<Pow>(ParType::DOUBLE, 2));
  Inverse->addNode(Denominator);
  Denominator->addNode(createLeaf(width));

  return tr;
}

struct PullInfo {
  double Mean;
  double MeanError;
  double Width;
  double WidthError;
};

PullInfo calculatePull(const std::vector<std::pair<double, double>> &FitValues,
                       double TrueValue) {

  std::vector<double> PullDist;

  for (auto const x : FitValues) {
    PullDist.push_back((x.first - TrueValue) / x.second);
  }

  PullInfo PInfo;
  // determine mean and sigma of the pull distribution
  double tempsum(0.0);
  for (auto const x : PullDist) {
    tempsum += x;
  }
  unsigned int N(PullDist.size());
  PInfo.Mean = tempsum / N;

  tempsum = 0.0;
  for (auto const x : PullDist) {
    tempsum += std::pow(x - PInfo.Mean, 2);
  }
  PInfo.Width = std::sqrt(tempsum / (N - 1));

  PInfo.MeanError = PInfo.Width / std::sqrt(N);
  PInfo.WidthError = PInfo.Width / std::sqrt(2.0 * N); // taken from ROOT

  return PInfo;
}

BOOST_AUTO_TEST_SUITE(Estimator_MinLogLHEstimatorTest)

BOOST_AUTO_TEST_CASE(MinLogLHEstimator_GaussianModelFitTest) {
  ComPWA::Logging log("INFO", "output.log");
  // the reference normal distribution of mean and sigma values
  double mean(3.0);
  double sigma(0.1);

  // IMPORTANT: set the starting values of the fit not to close to the edge of
  // the domain
  std::pair<double, double> domain_range(mean - 10.0 * sigma,
                                         mean + 10.0 * sigma);

  ComPWA::Data::DataSet PhspSample;
  PhspSample.Data.insert(std::make_pair("x", std::vector<double>()));
  std::mt19937 mt_gen(123456);

  std::uniform_real_distribution<double> distribution(domain_range.first,
                                                      domain_range.second);
  // a sample size of 20000 is necessary to have a good normalization
  // and precise fit parameters
  for (unsigned int i = 0; i < 20000; ++i) {
    PhspSample.Data["x"].push_back(distribution(mt_gen));
    PhspSample.Weights.push_back(1.0);
  }

  auto Gauss = Gaussian(mean, sigma);

  std::uniform_real_distribution<double> start_distribution(0.7, 1.3);
  std::normal_distribution<double> normal_distribution(mean, sigma);

  ComPWA::Data::DataSet DataSample;
  DataSample.Data.insert(std::make_pair("x", std::vector<double>()));
  DataSample.Weights = std::vector<double>(500, 1.0);
  std::chrono::milliseconds MeanFittime(0);
  std::chrono::milliseconds MeanFittimeFT(0);
  std::vector<std::pair<double, double>> MeanFitValues;
  std::vector<std::pair<double, double>> WidthFitValues;
  std::vector<std::pair<double, double>> MeanFitValuesFT;
  std::vector<std::pair<double, double>> WidthFitValuesFT;

  unsigned int NumSamples(50);
  for (unsigned i = 0; i < NumSamples; ++i) {
    DataSample.Data["x"].clear();
    for (unsigned int j = 0; j < 500; ++j) {
      DataSample.Data["x"].push_back(normal_distribution(mt_gen));
    }

    double startmean(start_distribution(mt_gen) * mean);
    double startsigma(start_distribution(mt_gen) * sigma);

    // this estimator is deprecated and will be removed soon
    ComPWA::FitParameterList InitialParameters;
    InitialParameters.push_back(
        ComPWA::FitParameter<double>("mean", startmean, false));
    InitialParameters.push_back(
        ComPWA::FitParameter<double>("width", startsigma, false));
    InitialParameters.push_back(
        ComPWA::FitParameter<double>("strength", 1.0, true));

    auto minLogLH = ComPWA::Estimator::MinLogLH(Gauss, DataSample, PhspSample);
    auto minuitif = ComPWA::Optimizer::Minuit2::MinuitIF();

    std::chrono::steady_clock::time_point StartTime =
        std::chrono::steady_clock::now();
    // STARTING MINIMIZATION
    auto result = minuitif.optimize(minLogLH, InitialParameters);
    std::chrono::steady_clock::time_point EndTime =
        std::chrono::steady_clock::now();

    MeanFittime += std::chrono::duration_cast<std::chrono::milliseconds>(
        EndTime - StartTime);
    MeanFitValues.push_back(
        std::make_pair(result.FinalParameters[0].Value,
                       result.FinalParameters[0].Error.first));
    WidthFitValues.push_back(
        std::make_pair(result.FinalParameters[1].Value,
                       result.FinalParameters[1].Error.first));

    BOOST_CHECK(std::abs(result.FinalParameters[0].Value - mean) <
                5.0 * result.FinalParameters[0].Error.first);
    BOOST_CHECK(std::abs(result.FinalParameters[1].Value - sigma) <
                5.0 * result.FinalParameters[1].Error.first);

    LOG(INFO) << "Now fitting using the function tree feature!";

    auto Mean =
        std::make_shared<ComPWA::FunctionTree::FitParameter>("Mean", startmean);
    Mean->fixParameter(false);
    auto Width = std::make_shared<ComPWA::FunctionTree::FitParameter>(
        "Width", startsigma);
    Width->fixParameter(false);
    auto Strength =
        std::make_shared<ComPWA::FunctionTree::FitParameter>("Strength", 1.0);
    ComPWA::FunctionTree::ParameterList Parameters;
    Strength->fixParameter(false);
    Parameters.addParameter(Mean);
    Parameters.addParameter(Width);
    Parameters.addParameter(Strength);

    ComPWA::FunctionTree::ParameterList DataList;
    DataList.addValue(
        std::make_shared<ComPWA::FunctionTree::Value<std::vector<double>>>(
            "x", std::vector<double>()));

    auto GaussFT = createFunctionTree(Mean, Width, Strength, DataList);

    auto intens = ComPWA::FunctionTree::FunctionTreeIntensity(
        GaussFT, Parameters, DataList);

    auto FTMinLogLH =
        ComPWA::Estimator::MinLogLH(intens, DataSample, PhspSample);

    minuitif = ComPWA::Optimizer::Minuit2::MinuitIF();

    std::chrono::steady_clock::time_point StartTimeFT =
        std::chrono::steady_clock::now();
    // STARTING MINIMIZATION
    auto resultft = minuitif.optimize(FTMinLogLH, InitialParameters);

    std::chrono::steady_clock::time_point EndTimeFT =
        std::chrono::steady_clock::now();

    MeanFittimeFT += std::chrono::duration_cast<std::chrono::milliseconds>(
        EndTimeFT - StartTimeFT);
    MeanFitValuesFT.push_back(
        std::make_pair(resultft.FinalParameters[0].Value,
                       resultft.FinalParameters[0].Error.first));
    WidthFitValuesFT.push_back(
        std::make_pair(resultft.FinalParameters[1].Value,
                       resultft.FinalParameters[1].Error.first));

    BOOST_CHECK(std::abs(resultft.FinalParameters[0].Value - mean) <
                5.0 * resultft.FinalParameters[0].Error.first);
    BOOST_CHECK(std::abs(resultft.FinalParameters[1].Value - sigma) <
                5.0 * resultft.FinalParameters[1].Error.first);
  }

  auto pm = calculatePull(MeanFitValues, mean);
  auto pw = calculatePull(WidthFitValues, sigma);
  auto pmft = calculatePull(MeanFitValuesFT, mean);
  auto pwft = calculatePull(WidthFitValuesFT, sigma);

  LOG(INFO) << "Mean pull: " << pm.Mean << "+-" << pm.MeanError << " | "
            << pm.Width << "+-" << pm.WidthError;
  LOG(INFO) << "Width pull: " << pw.Mean << "+-" << pw.MeanError << " | "
            << pw.Width << "+-" << pw.WidthError;

  LOG(INFO) << "Mean pull (FT): " << pmft.Mean << "+-" << pmft.MeanError
            << " | " << pmft.Width << "+-" << pmft.WidthError;
  LOG(INFO) << "Width pull (FT): " << pwft.Mean << "+-" << pwft.MeanError
            << " | " << pwft.Width << "+-" << pwft.WidthError;

  LOG(INFO) << "Mean fit runtime (w/o function tree): "
            << MeanFittime.count() / NumSamples << " ms";
  LOG(INFO) << "Mean fit runtime (with function tree): "
            << MeanFittimeFT.count() / NumSamples << " ms";

  // a pull should have mean == 0 and sigma == 1 (3 sigma error interval check)
  BOOST_CHECK(std::abs(pm.Mean) < 3.0 * pm.MeanError);
  BOOST_CHECK(std::abs(pm.Width - 1.0) < 3.0 * pm.WidthError);
  BOOST_CHECK(std::abs(pw.Mean) < 3.0 * pw.MeanError);
  BOOST_CHECK(std::abs(pw.Width - 1.0) < 3.0 * pw.WidthError);

  BOOST_CHECK(std::abs(pmft.Mean) < 3.0 * pmft.MeanError);
  BOOST_CHECK(std::abs(pmft.Width - 1.0) < 3.0 * pmft.WidthError);
  BOOST_CHECK(std::abs(pwft.Mean) < 3.0 * pwft.MeanError);
  BOOST_CHECK(std::abs(pwft.Width - 1.0) < 3.0 * pwft.WidthError);
};

BOOST_AUTO_TEST_CASE(MinLogLHEstimator_GaussianModelEventWeightTest) {
  ComPWA::Logging log("INFO", "output.log");
  // the reference normal distribution of mean and sigma values
  double mean(3.0);
  double sigma(0.1);

  // IMPORTANT: set the starting values of the fit not to close to the edge of
  // the domain
  std::pair<double, double> domain_range(mean - 10.0 * sigma,
                                         mean + 10.0 * sigma);

  ComPWA::Data::DataSet PhspSample;
  PhspSample.Data.insert(std::make_pair("x", std::vector<double>()));
  std::mt19937 mt_gen(123456);

  std::uniform_real_distribution<double> distribution(domain_range.first,
                                                      domain_range.second);
  // a sample size of 20000 is necessary to have a good normalization
  // and precise fit parameters
  for (unsigned int i = 0; i < 20000; ++i) {
    PhspSample.Data["x"].push_back(distribution(mt_gen));
    PhspSample.Weights.push_back(1.0);
  }

  auto Gauss = Gaussian(mean, sigma);

  // the integral needs to be calculated to normalize the Gaussian to set
  // appropriate weights for the data events
  // it is important to give the phase space volume, here the range of the
  // uniform number generation
  double integral = ComPWA::Tools::integrate(
      Gauss, PhspSample, domain_range.second - domain_range.first);
  LOG(INFO) << "Calculated integral: " << integral;

  ComPWA::FitParameterList InitialParameters;
  InitialParameters.push_back(ComPWA::FitParameter<double>("mean", 1.0, false));
  InitialParameters.push_back(
      ComPWA::FitParameter<double>("width", 1.0, false));
  InitialParameters.push_back(
      ComPWA::FitParameter<double>("strength", 1.0 / integral, true));

  std::uniform_real_distribution<double> start_distribution(0.7 * mean,
                                                            1.3 * mean);
  std::uniform_real_distribution<double> data_distribution(mean - 5.0 * sigma,
                                                           mean + 5.0 * sigma);
  ComPWA::Data::DataSet DataSample;
  DataSample.Data.insert(std::make_pair("x", std::vector<double>()));
  std::chrono::milliseconds MeanFittime(0);
  std::chrono::milliseconds MeanFittimeFT(0);
  std::vector<std::pair<double, double>> MeanFitValues;
  std::vector<std::pair<double, double>> WidthFitValues;
  std::vector<std::pair<double, double>> MeanFitValuesFT;
  std::vector<std::pair<double, double>> WidthFitValuesFT;

  unsigned int NumSamples(50);
  for (unsigned i = 0; i < NumSamples; ++i) {
    // reset the parameters, so that the weights are determined correctly
    Gauss.updateParametersFrom({mean, sigma, 1.0 / integral});

    DataSample.Data["x"].clear();
    DataSample.Weights.clear();
    unsigned int SampleSize(500);
    for (unsigned int j = 0; j < SampleSize; ++j) {
      DataSample.Data["x"].push_back(data_distribution(mt_gen));
    }
    auto TempIntensities = Gauss.evaluate(DataSample.Data);
    double WeightSum(
        std::accumulate(TempIntensities.begin(), TempIntensities.end(), 0.0));
    // reweight whole data sample so that weight sum = sample size
    // if this is not done, then the error estimates of the fit are
    // incorrect
    double RescaleFactor(1.0 * SampleSize / WeightSum);
    for (auto &x : TempIntensities) {
      x *= RescaleFactor;
    }
    DataSample.Weights = TempIntensities;

    double startmean(start_distribution(mt_gen));
    double startsigma(start_distribution(mt_gen) / mean * sigma);

    InitialParameters.at(0).Value = startmean;
    InitialParameters.at(1).Value = startsigma;

    LOG(INFO) << "Using start parameters " << startmean << " and "
              << startsigma;
    ComPWA::Estimator::MinLogLH minLogLH(Gauss, DataSample, PhspSample);

    auto minuitif = ComPWA::Optimizer::Minuit2::MinuitIF();

    std::chrono::steady_clock::time_point StartTime =
        std::chrono::steady_clock::now();
    // STARTING MINIMIZATION
    auto result = minuitif.optimize(minLogLH, InitialParameters);
    std::chrono::steady_clock::time_point EndTime =
        std::chrono::steady_clock::now();

    MeanFittime += std::chrono::duration_cast<std::chrono::milliseconds>(
        EndTime - StartTime);
    MeanFitValues.push_back(
        std::make_pair(result.FinalParameters[0].Value,
                       result.FinalParameters[0].Error.first));
    WidthFitValues.push_back(
        std::make_pair(result.FinalParameters[1].Value,
                       result.FinalParameters[1].Error.first));

    BOOST_CHECK(std::abs(result.FinalParameters[0].Value - mean) <
                5.0 * result.FinalParameters[0].Error.first);
    BOOST_CHECK(std::abs(std::abs(result.FinalParameters[1].Value) - sigma) <
                5.0 * result.FinalParameters[1].Error.first);

    LOG(INFO) << "Now fitting using the function tree feature!";

    auto Mean =
        std::make_shared<ComPWA::FunctionTree::FitParameter>("Mean", mean);
    Mean->fixParameter(false);
    auto Width =
        std::make_shared<ComPWA::FunctionTree::FitParameter>("Width", sigma);
    Width->fixParameter(false);
    auto Strength =
        std::make_shared<ComPWA::FunctionTree::FitParameter>("Strength", 1.0);
    ComPWA::FunctionTree::ParameterList Parameters;
    Strength->fixParameter(false);
    Parameters.addParameter(Mean);
    Parameters.addParameter(Width);
    Parameters.addParameter(Strength);

    ComPWA::FunctionTree::ParameterList DataList;
    DataList.addValue(
        std::make_shared<ComPWA::FunctionTree::Value<std::vector<double>>>(
            "x", std::vector<double>()));

    auto GaussFT = createFunctionTree(Mean, Width, Strength, DataList);

    auto intens = ComPWA::FunctionTree::FunctionTreeIntensity(
        GaussFT, Parameters, DataList);

    auto FTMinLogLH =
        ComPWA::Estimator::MinLogLH(intens, DataSample, PhspSample);

    std::chrono::steady_clock::time_point StartTimeFT =
        std::chrono::steady_clock::now();
    // STARTING MINIMIZATION
    auto resultft = minuitif.optimize(FTMinLogLH, InitialParameters);
    std::chrono::steady_clock::time_point EndTimeFT =
        std::chrono::steady_clock::now();

    MeanFittimeFT += std::chrono::duration_cast<std::chrono::milliseconds>(
        EndTimeFT - StartTimeFT);
    MeanFitValuesFT.push_back(
        std::make_pair(resultft.FinalParameters[0].Value,
                       resultft.FinalParameters[0].Error.first));
    WidthFitValuesFT.push_back(
        std::make_pair(resultft.FinalParameters[1].Value,
                       resultft.FinalParameters[1].Error.first));

    BOOST_CHECK(std::abs(resultft.FinalParameters[0].Value - mean) <
                5.0 * resultft.FinalParameters[0].Error.first);
    BOOST_CHECK(std::abs(std::abs(resultft.FinalParameters[1].Value) - sigma) <
                5.0 * resultft.FinalParameters[1].Error.first);
  }

  auto pm = calculatePull(MeanFitValues, mean);
  auto pw = calculatePull(WidthFitValues, sigma);
  auto pmft = calculatePull(MeanFitValuesFT, mean);
  auto pwft = calculatePull(WidthFitValuesFT, sigma);

  LOG(INFO) << "Mean pull: " << pm.Mean << "+-" << pm.MeanError << " | "
            << pm.Width << "+-" << pm.WidthError;
  LOG(INFO) << "Width pull: " << pw.Mean << "+-" << pw.MeanError << " | "
            << pw.Width << "+-" << pw.WidthError;

  LOG(INFO) << "Mean pull (FT): " << pmft.Mean << "+-" << pmft.MeanError
            << " | " << pmft.Width << "+-" << pmft.WidthError;
  LOG(INFO) << "Width pull (FT): " << pwft.Mean << "+-" << pwft.MeanError
            << " | " << pwft.Width << "+-" << pwft.WidthError;

  LOG(INFO) << "Mean fit runtime (w/o function tree): "
            << MeanFittime.count() / NumSamples << " ms";
  LOG(INFO) << "Mean fit runtime (with function tree): "
            << MeanFittimeFT.count() / NumSamples << " ms";

  // a pull should have mean == 0 and sigma == 1 (3 sigma error interval check)
  BOOST_CHECK(std::abs(pm.Mean) < 3.0 * pm.MeanError);
  BOOST_CHECK(std::abs(pm.Width - 1.0) < 3.0 * pm.WidthError);
  BOOST_CHECK(std::abs(pw.Mean) < 3.0 * pw.MeanError);
  BOOST_CHECK(std::abs(pw.Width - 1.0) < 3.0 * pw.WidthError);

  BOOST_CHECK(std::abs(pmft.Mean) < 3.0 * pmft.MeanError);
  BOOST_CHECK(std::abs(pmft.Width - 1.0) < 3.0 * pmft.WidthError);
  BOOST_CHECK(std::abs(pwft.Mean) < 3.0 * pwft.MeanError);
  BOOST_CHECK(std::abs(pwft.Width - 1.0) < 3.0 * pwft.WidthError);
};

BOOST_AUTO_TEST_SUITE_END()
