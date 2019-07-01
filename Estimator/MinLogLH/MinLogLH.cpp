// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "MinLogLH.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Intensity.hpp"
#include "Core/Kinematics.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Particle.hpp"
#include "Data/DataSet.hpp"

namespace ComPWA {
namespace Estimator {

MinLogLH::MinLogLH(std::shared_ptr<ComPWA::Intensity> intensity,
                   const std::vector<ComPWA::DataPoint> &datapoints,
                   const std::vector<ComPWA::DataPoint> &phsppoints)
    : Intensity(intensity), DataPoints(datapoints), PhspDataPoints(phsppoints) {

  LOG(INFO) << "MinLogLH::Init() |  Size of data sample = "
            << DataPoints.size();
}

double MinLogLH::evaluate() {
  double lh(0.0);

  double Norm(0.0);
  if (0 < PhspDataPoints.size()) {
    double PhspIntegral(0.0);
    double WeightSum(0.0);
    for (const auto &dp : PhspDataPoints) {
      PhspIntegral += Intensity->evaluate(dp.KinematicVariableList) * dp.Weight;
      WeightSum += dp.Weight;
    }

    Norm = (std::log(PhspIntegral / WeightSum) * DataPoints.size());
  }
  // calulate data log sum
  double LogSum(0.0);
  for (const auto &dp : DataPoints) {
    LogSum +=
        dp.Weight * std::log(Intensity->evaluate(dp.KinematicVariableList));
  }
  lh = Norm - LogSum;

  return lh;
}

std::shared_ptr<FunctionTree> createMinLogLHEstimatorFunctionTree(
    std::shared_ptr<ComPWA::OldIntensity> Intensity,
    std::shared_ptr<ComPWA::Data::DataSet> DataSample,
    std::shared_ptr<ComPWA::Data::DataSet> PhspDataSample) {
  LOG(DEBUG)
      << "createMinLogLHEstimatorFunctionTree(): constructing FunctionTree!";

  if (!DataSample) {
    LOG(ERROR)
        << "createMinLogLHEstimatorFunctionTree(): data sample is not set!";
    return {};
  }

  auto DataSampleList = DataSample->getParameterList();

  if (0 == DataSampleList.mDoubleValues().size()) {
    LOG(ERROR) << "createMinLogLHEstimatorFunctionTree(): Data sample is "
                  "empty! Please supply some data.";
    return std::shared_ptr<FunctionTree>(nullptr);
  }
  size_t SampleSize = DataSampleList.mDoubleValue(0)->values().size();

  std::shared_ptr<Value<std::vector<double>>> weights;
  try {
    weights = findMDoubleValue("Weight", DataSampleList);
  } catch (const Exception &e) {
  }

  std::shared_ptr<FunctionTree> EvaluationTree =
      std::make_shared<FunctionTree>("LH", std::make_shared<Value<double>>(),
                                     std::make_shared<AddAll>(ParType::DOUBLE));

  auto dataTree = std::make_shared<FunctionTree>(
      "DataEvaluation", std::make_shared<Value<double>>(),
      std::make_shared<MultAll>(ParType::DOUBLE));
  dataTree->createLeaf("minusOne", -1, "DataEvaluation");
  dataTree->createNode("Sum",
                       std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                       "DataEvaluation");
  dataTree->createNode("WeightedLogIntensities",
                       std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)),
                       "Sum");
  if (weights)
    dataTree->createLeaf("EventWeight", weights, "WeightedLogIntensities");
  dataTree->createNode("Log",
                       std::shared_ptr<Strategy>(new LogOf(ParType::MDOUBLE)),
                       "WeightedLogIntensities");
  dataTree->insertTree(Intensity->createFunctionTree(DataSampleList, ""),
                       "Log");

  EvaluationTree->insertTree(dataTree, "LH");

  // if there is a phasespace sample then do the normalization
  if (PhspDataSample) {
    auto PhspDataSampleList = PhspDataSample->getParameterList();
    if (0 < PhspDataSampleList.mDoubleValues().size()) {
      double PhspWeightSum(PhspDataSampleList.mDoubleValue(0)->values().size());

      std::shared_ptr<Value<std::vector<double>>> phspweights;
      try {
        phspweights = findMDoubleValue("Weight", PhspDataSampleList);
        PhspWeightSum = std::accumulate(phspweights->values().begin(),
                                        phspweights->values().end(), 0.0);
      } catch (const Exception &e) {
      }

      auto normTree = std::make_shared<FunctionTree>(
          "Normalization(intensity)", std::make_shared<Value<double>>(),
          std::make_shared<MultAll>(ParType::DOUBLE));
      normTree->createLeaf("N", SampleSize, "Normalization(intensity)");
      normTree->createNode(
          "Log", std::shared_ptr<Strategy>(new LogOf(ParType::DOUBLE)),
          "Normalization(intensity)");
      normTree->createNode(
          "Integral", std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)),
          "Log");
      // normTree->createLeaf("PhspVolume", PhspVolume, "Integral");
      normTree->createLeaf("InverseSampleWeights", 1.0 / PhspWeightSum,
                           "Integral");
      normTree->createNode(
          "Sum", std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
          "Integral");
      normTree->createNode(
          "WeightedIntensities",
          std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)), "Sum");
      if (phspweights)
        normTree->createLeaf("EventWeight", phspweights, "WeightedIntensities");
      normTree->insertTree(
          Intensity->createFunctionTree(PhspDataSampleList, "phsp"),
          "WeightedIntensities");

      EvaluationTree->insertTree(normTree, "LH");
    }
  } else {
    LOG(INFO) << "createMinLogLHEstimatorFunctionTree(): phsp sample is empty! "
                 "Skipping normalization and assuming intensity is normalized!";
  }
  LOG(DEBUG) << "createMinLogLHEstimatorFunctionTree(): construction of LH "
                "tree finished! Performing checks ...";
  EvaluationTree->parameter();
  if (!EvaluationTree->sanityCheck()) {
    throw std::runtime_error(
        "createMinLogLHEstimatorFunctionTree(): tree has structural "
        "problems. Sanity check not passed!");
  }
  LOG(DEBUG) << "createMinLogLHEstimatorFunctionTree(): finished!";
  return EvaluationTree;
}

std::shared_ptr<FunctionTreeEstimatorWrapper>
createMinLogLHFunctionTreeEstimator(
    std::shared_ptr<ComPWA::OldIntensity> Intensity,
    std::shared_ptr<ComPWA::Data::DataSet> DataSample,
    std::shared_ptr<ComPWA::Data::DataSet> PhspDataSample) {

  auto ft = createMinLogLHEstimatorFunctionTree(Intensity, DataSample,
                                                PhspDataSample);

  ParameterList params;
  Intensity->addUniqueParametersTo(params);

  return std::make_shared<FunctionTreeEstimatorWrapper>(ft, params);
}

} // namespace Estimator
} // namespace ComPWA
