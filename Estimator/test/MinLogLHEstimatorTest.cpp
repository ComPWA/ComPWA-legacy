// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Estimator_MinLogLHEstimatorTest

#include <chrono>
#include <cmath>
#include <random>

#include <boost/test/unit_test.hpp>

#include "Core/Event.hpp"
#include "Core/FunctionTreeIntensityWrapper.hpp"
#include "Core/Intensity.hpp"
#include "Core/ParameterList.hpp"
#include "Data/DataSet.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"
#include "Tools/Integration.hpp"

BOOST_AUTO_TEST_SUITE(Estimator_MinLogLHEstimatorTest)

class Gaussian : public ComPWA::OldIntensity {
  std::shared_ptr<ComPWA::FitParameter> Mean;
  std::shared_ptr<ComPWA::FitParameter> Width;

  std::shared_ptr<ComPWA::FitParameter> Strength;

public:
  Gaussian(double mean, double width) {
    Mean = std::make_shared<ComPWA::FitParameter>("Mean", mean);
    Width = std::make_shared<ComPWA::FitParameter>("Width", width);
    Strength = std::make_shared<ComPWA::FitParameter>("Strength", 1.0);
  }

  double evaluate(const ComPWA::DataPoint &point) const {
    return Strength->value() *
           std::exp(
               -0.5 *
               std::pow(point.KinematicVariableList[0] - Mean->value(), 2) /
               std::pow(Width->value(), 2));
  }

  std::shared_ptr<ComPWA::FunctionTree>
  createFunctionTree(const ComPWA::ParameterList &DataSample,
                     const std::string &suffix) const {

    size_t n = DataSample.mDoubleValue(0)->values().size();

    std::string Name("Gaussian");
    auto tr = std::make_shared<ComPWA::FunctionTree>(
        Name + suffix, ComPWA::MDouble("", n),
        std::shared_ptr<ComPWA::Strategy>(
            new ComPWA::MultAll(ComPWA::ParType::MDOUBLE)));
    tr->createLeaf("Strength", Strength, Name + suffix);
    tr->createNode("Exp", ComPWA::MDouble("", n),
                   std::make_shared<ComPWA::Exp>(ComPWA::ParType::MDOUBLE),
                   Name + suffix);
    tr->createNode("Exponent", ComPWA::MDouble("", n),
                   std::make_shared<ComPWA::MultAll>(ComPWA::ParType::MDOUBLE),
                   "Exp");
    tr->createLeaf("minusHalf", -0.5, "Exponent");

    tr->createNode("Nominator", ComPWA::MDouble("", n),
                   std::make_shared<ComPWA::Pow>(ComPWA::ParType::MDOUBLE, 2),
                   "Exponent");
    tr->createNode("Diff",
                   std::make_shared<ComPWA::AddAll>(ComPWA::ParType::MDOUBLE),
                   "Nominator");
    tr->createLeaf("x", DataSample.mDoubleValue(0), "Diff");
    tr->createNode("Negate",
                   std::shared_ptr<ComPWA::Parameter>(
                       new ComPWA::Value<double>("negMean")),
                   std::make_shared<ComPWA::MultAll>(ComPWA::ParType::DOUBLE),
                   "Diff");
    tr->createLeaf("minusOne", -1.0, "Negate");
    tr->createLeaf("Mean", Mean, "Negate");

    tr->createNode("Inverse",
                   std::shared_ptr<ComPWA::Parameter>(
                       new ComPWA::Value<double>("invdenom")),
                   std::make_shared<ComPWA::Inverse>(ComPWA::ParType::DOUBLE),
                   "Exponent");
    tr->createNode(
        "Denominator",
        std::shared_ptr<ComPWA::Parameter>(new ComPWA::Value<double>("denom")),
        std::make_shared<ComPWA::Pow>(ComPWA::ParType::DOUBLE, 2), "Inverse");
    tr->createLeaf("Width", Width, "Denominator");

    return tr;
  }

  void addUniqueParametersTo(ComPWA::ParameterList &list) {
    Strength = list.addUniqueParameter(Strength);
    Mean = list.addUniqueParameter(Mean);
    Width = list.addUniqueParameter(Width);
  }

  void addFitParametersTo(std::vector<double> &list) {
    list.push_back(Strength->value());
    list.push_back(Mean->value());
    list.push_back(Width->value());
  }

  void updateParametersFrom(const ComPWA::ParameterList &list) {
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

  std::vector<ComPWA::DataPoint> PhspDataPoints;
  std::mt19937 mt_gen(123456);

  std::uniform_real_distribution<double> distribution(domain_range.first,
                                                      domain_range.second);
  // a sample size of 20000 is necessary to have a good normalization
  // and precise fit parameters
  for (unsigned int i = 0; i < 20000; ++i) {
    ComPWA::DataPoint dp;
    dp.KinematicVariableList.push_back(distribution(mt_gen));
    PhspDataPoints.push_back(dp);
  }
  auto PhspSample = std::make_shared<ComPWA::Data::DataSet>(PhspDataPoints);

  std::shared_ptr<ComPWA::OldIntensity> Gauss(new Gaussian(mean, sigma));

  ComPWA::ParameterList FitParameters;
  Gauss->addUniqueParametersTo(FitParameters);

  auto MeanParameter = ComPWA::FindParameter("Mean", FitParameters);
  MeanParameter->fixParameter(false);
  auto WidthParameter = ComPWA::FindParameter("Width", FitParameters);
  WidthParameter->fixParameter(false);

  std::uniform_real_distribution<double> start_distribution(0.7, 1.3);
  std::normal_distribution<double> normal_distribution(mean, sigma);

  std::vector<ComPWA::DataPoint> DataPoints;
  std::chrono::milliseconds MeanFittime(0);
  std::chrono::milliseconds MeanFittimeFT(0);
  std::vector<std::pair<double, double>> MeanFitValues;
  std::vector<std::pair<double, double>> WidthFitValues;
  std::vector<std::pair<double, double>> MeanFitValuesFT;
  std::vector<std::pair<double, double>> WidthFitValuesFT;

  unsigned int NumSamples(50);
  for (unsigned i = 0; i < NumSamples; ++i) {
    DataPoints.clear();
    for (unsigned int j = 0; j < 500; ++j) {
      ComPWA::DataPoint dp;
      dp.KinematicVariableList.push_back(normal_distribution(mt_gen));
      DataPoints.push_back(dp);
    }
    auto DataSample = std::make_shared<ComPWA::Data::DataSet>(DataPoints);

    double startmean(start_distribution(mt_gen) * mean);
    double startsigma(start_distribution(mt_gen) * sigma);

    MeanParameter->setValue(startmean);
    WidthParameter->setValue(startsigma);

    auto intens = std::make_shared<ComPWA::FunctionTreeIntensityWrapper>(
        Gauss, 1, "gauss");

    auto minLogLH = std::make_shared<ComPWA::Estimator::MinLogLH>(
        intens, DataSample->getParameterList(), PhspSample->getParameterList());

    auto minuitif =
        new ComPWA::Optimizer::Minuit2::MinuitIF(minLogLH, FitParameters);
    minuitif->setUseHesse(true);

    std::chrono::steady_clock::time_point StartTime =
        std::chrono::steady_clock::now();
    // STARTING MINIMIZATION
    auto result =
        std::dynamic_pointer_cast<ComPWA::Optimizer::Minuit2::MinuitResult>(
            minuitif->exec(FitParameters));
    std::chrono::steady_clock::time_point EndTime =
        std::chrono::steady_clock::now();

    MeanFittime += std::chrono::duration_cast<std::chrono::milliseconds>(
        EndTime - StartTime);
    MeanFitValues.push_back(
        std::make_pair(MeanParameter->value(), MeanParameter->error().first));
    WidthFitValues.push_back(
        std::make_pair(WidthParameter->value(), WidthParameter->error().first));

    BOOST_CHECK(std::abs(MeanParameter->value() - mean) <
                5.0 * MeanParameter->error().first);
    BOOST_CHECK(std::abs(WidthParameter->value() - sigma) <
                5.0 * WidthParameter->error().first);

    LOG(INFO) << "Now fitting using the function tree feature!";

    MeanParameter->setValue(startmean);
    WidthParameter->setValue(startsigma);
    MeanParameter->setError(0.0);
    WidthParameter->setError(0.0);

    auto FTMinLogLH = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
        Gauss, DataSample, PhspSample);

    minuitif =
        new ComPWA::Optimizer::Minuit2::MinuitIF(FTMinLogLH, FitParameters);
    minuitif->setUseHesse(true);

    std::chrono::steady_clock::time_point StartTimeFT =
        std::chrono::steady_clock::now();
    // STARTING MINIMIZATION
    result =
        std::dynamic_pointer_cast<ComPWA::Optimizer::Minuit2::MinuitResult>(
            minuitif->exec(FitParameters));
    std::chrono::steady_clock::time_point EndTimeFT =
        std::chrono::steady_clock::now();

    MeanFittimeFT += std::chrono::duration_cast<std::chrono::milliseconds>(
        EndTimeFT - StartTimeFT);
    MeanFitValuesFT.push_back(
        std::make_pair(MeanParameter->value(), MeanParameter->error().first));
    WidthFitValuesFT.push_back(
        std::make_pair(WidthParameter->value(), WidthParameter->error().first));

    BOOST_CHECK(std::abs(MeanParameter->value() - mean) <
                5.0 * MeanParameter->error().first);
    BOOST_CHECK(std::abs(WidthParameter->value() - sigma) <
                5.0 * WidthParameter->error().first);
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

  std::vector<ComPWA::DataPoint> PhspDataPoints;
  std::mt19937 mt_gen(123456);

  std::uniform_real_distribution<double> distribution(domain_range.first,
                                                      domain_range.second);
  // a sample size of 20000 is necessary to have a good normalization
  // and precise fit parameters
  for (unsigned int i = 0; i < 20000; ++i) {
    ComPWA::DataPoint dp;
    dp.KinematicVariableList.push_back(distribution(mt_gen));
    PhspDataPoints.push_back(dp);
  }
  auto PhspSample = std::make_shared<ComPWA::Data::DataSet>(PhspDataPoints);

  std::shared_ptr<ComPWA::OldIntensity> Gauss(new Gaussian(mean, sigma));

  // the integral needs to be calculated to normalize the Gaussian to set
  // appropriate weights for the data events
  // it is important to give the phase space volume, here the range of the
  // uniform number generation
  double integral = ComPWA::Tools::integrate(
      Gauss, PhspSample, domain_range.second - domain_range.first);
  LOG(INFO) << "Calculated integral: " << integral;

  ComPWA::ParameterList FitParameters;
  Gauss->addUniqueParametersTo(FitParameters);

  auto StrengthParameter = ComPWA::FindParameter("Strength", FitParameters);
  StrengthParameter->fixParameter(false);
  StrengthParameter->setValue(1.0 / integral);
  StrengthParameter->fixParameter(true);

  auto MeanParameter = ComPWA::FindParameter("Mean", FitParameters);
  MeanParameter->fixParameter(false);
  auto WidthParameter = ComPWA::FindParameter("Width", FitParameters);
  WidthParameter->fixParameter(false);

  std::uniform_real_distribution<double> start_distribution(0.7 * mean,
                                                            1.3 * mean);
  std::uniform_real_distribution<double> data_distribution(mean - 5.0 * sigma,
                                                           mean + 5.0 * sigma);
  std::vector<ComPWA::DataPoint> DataPoints;
  std::chrono::milliseconds MeanFittime(0);
  std::chrono::milliseconds MeanFittimeFT(0);
  std::vector<std::pair<double, double>> MeanFitValues;
  std::vector<std::pair<double, double>> WidthFitValues;
  std::vector<std::pair<double, double>> MeanFitValuesFT;
  std::vector<std::pair<double, double>> WidthFitValuesFT;

  unsigned int NumSamples(50);
  for (unsigned i = 0; i < NumSamples; ++i) {
    MeanParameter->setValue(mean);
    WidthParameter->setValue(sigma);

    DataPoints.clear();
    double WeightSum(0.0);
    unsigned int SampleSize(2000);
    for (unsigned int j = 0; j < SampleSize; ++j) {
      ComPWA::DataPoint dp;
      dp.KinematicVariableList.push_back(data_distribution(mt_gen));
      dp.Weight = Gauss->evaluate(dp);
      WeightSum += dp.Weight;
      DataPoints.push_back(dp);
    }
    // reweight whole data sample so that weight sum = sample size
    // if this is not done, then the error estimates of the fit are incorrect
    double RescaleFactor(1.0 * SampleSize / WeightSum);
    for (auto &x : DataPoints) {
      x.Weight = x.Weight * RescaleFactor;
    }

    auto DataSample = std::make_shared<ComPWA::Data::DataSet>(DataPoints);

    double startmean(start_distribution(mt_gen));
    double startsigma(start_distribution(mt_gen) / mean * sigma);

    MeanParameter->setValue(startmean);
    WidthParameter->setValue(startsigma);

    LOG(INFO) << "Using start parameters " << startmean << " and "
              << startsigma;
    auto intens = std::make_shared<ComPWA::FunctionTreeIntensityWrapper>(
        Gauss, 1, "gauss");
    std::shared_ptr<ComPWA::Estimator::MinLogLH> minLogLH(
        new ComPWA::Estimator::MinLogLH(intens, DataSample->getParameterList(),
                                        PhspSample->getParameterList()));

    auto minuitif =
        new ComPWA::Optimizer::Minuit2::MinuitIF(minLogLH, FitParameters);
    minuitif->setUseHesse(true);

    std::chrono::steady_clock::time_point StartTime =
        std::chrono::steady_clock::now();
    // STARTING MINIMIZATION
    auto result =
        std::dynamic_pointer_cast<ComPWA::Optimizer::Minuit2::MinuitResult>(
            minuitif->exec(FitParameters));
    std::chrono::steady_clock::time_point EndTime =
        std::chrono::steady_clock::now();

    MeanFittime += std::chrono::duration_cast<std::chrono::milliseconds>(
        EndTime - StartTime);
    MeanFitValues.push_back(
        std::make_pair(MeanParameter->value(), MeanParameter->error().first));
    WidthFitValues.push_back(
        std::make_pair(WidthParameter->value(), WidthParameter->error().first));

    BOOST_CHECK(std::abs(MeanParameter->value() - mean) <
                5.0 * MeanParameter->error().first);
    BOOST_CHECK(std::abs(std::abs(WidthParameter->value()) - sigma) <
                5.0 * WidthParameter->error().first);

    LOG(INFO) << "Now fitting using the function tree feature!";
    MeanParameter->setValue(startmean);
    WidthParameter->setValue(startsigma);
    MeanParameter->setError(0.0);
    WidthParameter->setError(0.0);

    auto FTMinLogLH = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
        Gauss, DataSample, PhspSample);

    minuitif =
        new ComPWA::Optimizer::Minuit2::MinuitIF(FTMinLogLH, FitParameters);
    minuitif->setUseHesse(true);

    std::chrono::steady_clock::time_point StartTimeFT =
        std::chrono::steady_clock::now();
    // STARTING MINIMIZATION
    result =
        std::dynamic_pointer_cast<ComPWA::Optimizer::Minuit2::MinuitResult>(
            minuitif->exec(FitParameters));
    std::chrono::steady_clock::time_point EndTimeFT =
        std::chrono::steady_clock::now();

    MeanFittimeFT += std::chrono::duration_cast<std::chrono::milliseconds>(
        EndTimeFT - StartTimeFT);
    MeanFitValuesFT.push_back(
        std::make_pair(MeanParameter->value(), MeanParameter->error().first));
    WidthFitValuesFT.push_back(
        std::make_pair(WidthParameter->value(), WidthParameter->error().first));

    BOOST_CHECK(std::abs(MeanParameter->value() - mean) <
                5.0 * MeanParameter->error().first);
    BOOST_CHECK(std::abs(std::abs(WidthParameter->value()) - sigma) <
                5.0 * WidthParameter->error().first);
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
