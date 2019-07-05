// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Estimator_MinLogLHEstimatorTest

#include <chrono>
#include <cmath>
#include <random>

#include <boost/test/unit_test.hpp>

#include "Core/FunctionTree/FunctionTreeIntensityWrapper.hpp"
#include "Core/FunctionTree/Intensity.hpp"
#include "Core/FunctionTree/ParameterList.hpp"
#include "Data/DataSet.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"
#include "Tools/Integration.hpp"

using namespace ComPWA::FunctionTree;

BOOST_AUTO_TEST_SUITE(Estimator_MinLogLHEstimatorTest)

class Gaussian : public ComPWA::FunctionTree::OldIntensity {
  std::shared_ptr<FitParameter> Mean;
  std::shared_ptr<FitParameter> Width;

  std::shared_ptr<FitParameter> Strength;

public:
  Gaussian(double mean, double width) {
    Mean = std::make_shared<FitParameter>("Mean", mean);
    Width = std::make_shared<FitParameter>("Width", width);
    Strength = std::make_shared<FitParameter>("Strength", 1.0);
  }

  double evaluate(const ComPWA::DataPoint &point) const {
    return Strength->value() *
           std::exp(
               -0.5 *
               std::pow(point.KinematicVariableList[0] - Mean->value(), 2) /
               std::pow(Width->value(), 2));
  }

  std::shared_ptr<ComPWA::FunctionTree::FunctionTree>
  createFunctionTree(const ParameterList &DataSample,
                     const std::string &suffix) const {

    size_t n = DataSample.mDoubleValue(0)->values().size();

    std::string Name("Gaussian");
    auto tr = std::make_shared<ComPWA::FunctionTree::FunctionTree>(
        Name + suffix, MDouble("", n),
        std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
    tr->createLeaf("Strength", Strength, Name + suffix);
    tr->createNode("Exp", MDouble("", n),
                   std::make_shared<Exp>(ParType::MDOUBLE), Name + suffix);
    tr->createNode("Exponent", MDouble("", n),
                   std::make_shared<MultAll>(ParType::MDOUBLE), "Exp");
    tr->createLeaf("minusHalf", -0.5, "Exponent");

    tr->createNode("Nominator", MDouble("", n),
                   std::make_shared<Pow>(ParType::MDOUBLE, 2), "Exponent");
    tr->createNode("Diff", std::make_shared<AddAll>(ParType::MDOUBLE),
                   "Nominator");
    tr->createLeaf("x", DataSample.mDoubleValue(0), "Diff");
    tr->createNode("Negate",
                   std::shared_ptr<Parameter>(new Value<double>("negMean")),
                   std::make_shared<MultAll>(ParType::DOUBLE), "Diff");
    tr->createLeaf("minusOne", -1.0, "Negate");
    tr->createLeaf("Mean", Mean, "Negate");

    tr->createNode("Inverse",
                   std::shared_ptr<Parameter>(new Value<double>("invdenom")),
                   std::make_shared<Inverse>(ParType::DOUBLE), "Exponent");
    tr->createNode("Denominator",
                   std::shared_ptr<Parameter>(new Value<double>("denom")),
                   std::make_shared<Pow>(ParType::DOUBLE, 2), "Inverse");
    tr->createLeaf("Width", Width, "Denominator");

    return tr;
  }

  void addUniqueParametersTo(ParameterList &list) {
    Strength = list.addUniqueParameter(Strength);
    Mean = list.addUniqueParameter(Mean);
    Width = list.addUniqueParameter(Width);
  }

  void addFitParametersTo(std::vector<double> &list) {
    list.push_back(Strength->value());
    list.push_back(Mean->value());
    list.push_back(Width->value());
  }

  void updateParametersFrom(const ParameterList &list) {
    auto p = FindParameter(Strength->name(), list);
    Strength->updateParameter(p);
    p = FindParameter(Mean->name(), list);
    Mean->updateParameter(p);
    p = FindParameter(Width->name(), list);
    Width->updateParameter(p);
  }
};

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

BOOST_AUTO_TEST_CASE(MinLogLHEstimator_GaussianModelFitTest) {
  ComPWA::Logging log("output.log", "INFO");
  // the reference normal distribution of mean and sigma values
  double mean(3.0);
  double sigma(0.1);

  // IMPORTANT: set the starting values of the fit not to close to the edge of
  // the domain
  std::pair<double, double> domain_range(mean - 10.0 * sigma,
                                         mean + 10.0 * sigma);

  ComPWA::Data::DataSet PhspSample;
  PhspSample.Data.push_back({});
  std::mt19937 mt_gen(123456);

  std::uniform_real_distribution<double> distribution(domain_range.first,
                                                      domain_range.second);
  // a sample size of 20000 is necessary to have a good normalization
  // and precise fit parameters
  for (unsigned int i = 0; i < 20000; ++i) {
    PhspSample.Data[0].push_back(distribution(mt_gen));
    PhspSample.Weights.push_back(1.0);
  }

  std::shared_ptr<OldIntensity> Gauss(new Gaussian(mean, sigma));

  ParameterList FitParameters;
  Gauss->addUniqueParametersTo(FitParameters);

  /*auto MeanParameter = FindParameter("Mean", FitParameters);
  MeanParameter->fixParameter(false);
  auto WidthParameter = FindParameter("Width", FitParameters);
  WidthParameter->fixParameter(false);*/

  std::uniform_real_distribution<double> start_distribution(0.7, 1.3);
  std::normal_distribution<double> normal_distribution(mean, sigma);

  ComPWA::Data::DataSet DataSample;
  DataSample.Data.push_back({});
  DataSample.Weights = std::vector<double>(500, 1.0);
  std::chrono::milliseconds MeanFittime(0);
  std::chrono::milliseconds MeanFittimeFT(0);
  std::vector<std::pair<double, double>> MeanFitValues;
  std::vector<std::pair<double, double>> WidthFitValues;
  std::vector<std::pair<double, double>> MeanFitValuesFT;
  std::vector<std::pair<double, double>> WidthFitValuesFT;

  unsigned int NumSamples(50);
  for (unsigned i = 0; i < NumSamples; ++i) {
    DataSample.Data[0].clear();
    for (unsigned int j = 0; j < 500; ++j) {
      DataSample.Data[0].push_back(normal_distribution(mt_gen));
    }

    double startmean(start_distribution(mt_gen) * mean);
    double startsigma(start_distribution(mt_gen) * sigma);

    auto intens =
        std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensityWrapper>(
            Gauss, 1, "gauss");

    // this estimator is deprecated and will be removed soon
    ComPWA::FitParameterList InitialParameters;
    InitialParameters.push_back(
        ComPWA::FitParameter<double>("strength", 1.0, true));
    InitialParameters.push_back(
        ComPWA::FitParameter<double>("mean", startmean, false));
    InitialParameters.push_back(
        ComPWA::FitParameter<double>("width", startsigma, false));

    auto minLogLH = ComPWA::Estimator::MinLogLH(intens, DataSample, PhspSample);
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
        std::make_pair(result.FinalParameters[1].Value,
                       result.FinalParameters[1].Error.first));
    WidthFitValues.push_back(
        std::make_pair(result.FinalParameters[2].Value,
                       result.FinalParameters[2].Error.first));

    BOOST_CHECK(std::abs(result.FinalParameters[1].Value - mean) <
                5.0 * result.FinalParameters[1].Error.first);
    BOOST_CHECK(std::abs(result.FinalParameters[2].Value - sigma) <
                5.0 * result.FinalParameters[2].Error.first);

    LOG(INFO) << "Now fitting using the function tree feature!";

    auto FTMinLogLH = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
        intens, DataSample, PhspSample);

    minuitif = ComPWA::Optimizer::Minuit2::MinuitIF();

    std::chrono::steady_clock::time_point StartTimeFT =
        std::chrono::steady_clock::now();
    // STARTING MINIMIZATION
    result = minuitif.optimize(std::get<0>(FTMinLogLH), InitialParameters);
    std::chrono::steady_clock::time_point EndTimeFT =
        std::chrono::steady_clock::now();

    MeanFittimeFT += std::chrono::duration_cast<std::chrono::milliseconds>(
        EndTimeFT - StartTimeFT);
    MeanFitValuesFT.push_back(
        std::make_pair(result.FinalParameters[1].Value,
                       result.FinalParameters[1].Error.first));
    WidthFitValuesFT.push_back(
        std::make_pair(result.FinalParameters[2].Value,
                       result.FinalParameters[2].Error.first));

    BOOST_CHECK(std::abs(result.FinalParameters[1].Value - mean) <
                5.0 * result.FinalParameters[1].Error.first);
    BOOST_CHECK(std::abs(result.FinalParameters[2].Value - sigma) <
                5.0 * result.FinalParameters[2].Error.first);
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
  ComPWA::Logging log("output.log", "INFO");
  // the reference normal distribution of mean and sigma values
  double mean(3.0);
  double sigma(0.1);

  // IMPORTANT: set the starting values of the fit not to close to the edge of
  // the domain
  std::pair<double, double> domain_range(mean - 10.0 * sigma,
                                         mean + 10.0 * sigma);

  ComPWA::Data::DataSet PhspSample;
  PhspSample.Data.push_back({});
  std::mt19937 mt_gen(123456);

  std::uniform_real_distribution<double> distribution(domain_range.first,
                                                      domain_range.second);
  // a sample size of 20000 is necessary to have a good normalization
  // and precise fit parameters
  for (unsigned int i = 0; i < 20000; ++i) {
    PhspSample.Data[0].push_back(distribution(mt_gen));
    PhspSample.Weights.push_back(1.0);
  }

  auto GaussOld = std::make_shared<Gaussian>(mean, sigma);

  auto Gauss =
      std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensityWrapper>(
          GaussOld, 1, "gauss");

  // the integral needs to be calculated to normalize the Gaussian to set
  // appropriate weights for the data events
  // it is important to give the phase space volume, here the range of the
  // uniform number generation
  double integral = ComPWA::Tools::integrate(
      Gauss, PhspSample, domain_range.second - domain_range.first);
  LOG(INFO) << "Calculated integral: " << integral;

  ComPWA::FitParameterList InitialParameters;
  InitialParameters.push_back(
      ComPWA::FitParameter<double>("strength", 1.0 / integral, true));
  InitialParameters.push_back(ComPWA::FitParameter<double>("mean", 1.0, false));
  InitialParameters.push_back(
      ComPWA::FitParameter<double>("width", 1.0, false));

  std::uniform_real_distribution<double> start_distribution(0.7 * mean,
                                                            1.3 * mean);
  std::uniform_real_distribution<double> data_distribution(mean - 5.0 * sigma,
                                                           mean + 5.0 * sigma);
  ComPWA::Data::DataSet DataSample;
  DataSample.Data.push_back({});
  std::chrono::milliseconds MeanFittime(0);
  std::chrono::milliseconds MeanFittimeFT(0);
  std::vector<std::pair<double, double>> MeanFitValues;
  std::vector<std::pair<double, double>> WidthFitValues;
  std::vector<std::pair<double, double>> MeanFitValuesFT;
  std::vector<std::pair<double, double>> WidthFitValuesFT;

  unsigned int NumSamples(50);
  for (unsigned i = 0; i < NumSamples; ++i) {
    // reset the parameters, so that the weights are determined correctly
    Gauss->updateParametersFrom({1.0 / integral, mean, sigma});

    DataSample.Data[0].clear();
    DataSample.Weights.clear();
    unsigned int SampleSize(500);
    for (unsigned int j = 0; j < SampleSize; ++j) {
      DataSample.Data[0].push_back(data_distribution(mt_gen));
    }
    auto TempIntensities = Gauss->evaluate(DataSample.Data);
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

    InitialParameters.at(1).Value = startmean;
    InitialParameters.at(2).Value = startsigma;

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
        std::make_pair(result.FinalParameters[1].Value,
                       result.FinalParameters[1].Error.first));
    WidthFitValues.push_back(
        std::make_pair(result.FinalParameters[2].Value,
                       result.FinalParameters[2].Error.first));

    BOOST_CHECK(std::abs(result.FinalParameters[1].Value - mean) <
                5.0 * result.FinalParameters[1].Error.first);
    BOOST_CHECK(std::abs(std::abs(result.FinalParameters[2].Value) - sigma) <
                5.0 * result.FinalParameters[2].Error.first);

    LOG(INFO) << "Now fitting using the function tree feature!";
    InitialParameters.at(1).Value = startmean;
    InitialParameters.at(2).Value = startsigma;

    auto FTMinLogLH = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
        Gauss, DataSample, PhspSample);

    std::chrono::steady_clock::time_point StartTimeFT =
        std::chrono::steady_clock::now();
    // STARTING MINIMIZATION
    result = minuitif.optimize(std::get<0>(FTMinLogLH), InitialParameters);
    std::chrono::steady_clock::time_point EndTimeFT =
        std::chrono::steady_clock::now();

    MeanFittimeFT += std::chrono::duration_cast<std::chrono::milliseconds>(
        EndTimeFT - StartTimeFT);
    MeanFitValuesFT.push_back(
        std::make_pair(result.FinalParameters[1].Value,
                       result.FinalParameters[1].Error.first));
    WidthFitValuesFT.push_back(
        std::make_pair(result.FinalParameters[2].Value,
                       result.FinalParameters[2].Error.first));

    BOOST_CHECK(std::abs(result.FinalParameters[1].Value - mean) <
                5.0 * result.FinalParameters[1].Error.first);
    BOOST_CHECK(std::abs(std::abs(result.FinalParameters[2].Value) - sigma) <
                5.0 * result.FinalParameters[2].Error.first);
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
