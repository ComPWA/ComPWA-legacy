// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "MinLogLH.hpp"
#include "Core/FunctionTree/ParameterList.hpp"
#include "Core/FunctionTree/TreeNode.hpp"
#include "Core/Kinematics.hpp"
#include "Core/FourMomentum.hpp"
#include "Data/DataSet.hpp"

namespace ComPWA {
namespace Estimator {

using namespace ComPWA::FunctionTree;

MinLogLH::MinLogLH(ComPWA::Intensity &intensity,
                   const Data::DataSet &datasample,
                   const Data::DataSet &phspdatasample)
    : Intensity(intensity), DataSample(datasample),
      PhspDataSample(phspdatasample) {

  LOG(INFO) << "MinLogLH::MinLogLH() |  Size of data sample = "
            << DataSample.Weights.size();
}

double MinLogLH::evaluate() noexcept {
  double lh(0.0);

  double Norm(0.0);
  if (0 < PhspDataSample.Weights.size()) {
    double PhspIntegral(0.0);
    double WeightSum(0.0);
    auto Intensities = Intensity.evaluate(PhspDataSample.Data);
    auto IntensIter = Intensities.begin();
    for (auto x = PhspDataSample.Weights.begin();
         x != PhspDataSample.Weights.end(); ++x) {
      PhspIntegral += *x * *IntensIter;
      WeightSum += *x;
      ++IntensIter;
    }
    Norm = (std::log(PhspIntegral / WeightSum) * DataSample.Weights.size());
  }
  // calculate data log sum
  double LogSum(0.0);
  auto Intensities = Intensity.evaluate(DataSample.Data);
  for (size_t i = 0; i < DataSample.Weights.size(); ++i) {
    LogSum += std::log(Intensities[i]) * DataSample.Weights[i];
  }
  lh = Norm - LogSum;

  return lh;
}

void MinLogLH::updateParametersFrom(const std::vector<double> &params) {
  Intensity.updateParametersFrom(params);
}

std::vector<ComPWA::Parameter> MinLogLH::getParameters() const {
  return Intensity.getParameters();
}

std::pair<ComPWA::FunctionTree::FunctionTreeEstimator, FitParameterList>
createMinLogLHFunctionTreeEstimator(
    ComPWA::FunctionTree::FunctionTreeIntensity &Intensity,
    const ComPWA::Data::DataSet &DataSample) {
  using namespace ComPWA::FunctionTree;
  LOG(DEBUG)
      << "createMinLogLHEstimatorFunctionTree(): constructing FunctionTree!";

  if (0 == DataSample.Weights.size()) {
    LOG(ERROR) << "createMinLogLHEstimatorFunctionTree(): Data sample is "
                  "empty! Please supply some data.";
    return std::make_pair(FunctionTreeEstimator(nullptr, ParameterList()),
                          FitParameterList());
  }

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> DataIntensityFunctionTree;
  ComPWA::FunctionTree::ParameterList Parameters;
  std::tie(DataIntensityFunctionTree, Parameters) =
      Intensity.bind(DataSample.Data);

  FitParameterList FitParList =
      ComPWA::FunctionTree::createFitParameterList(Parameters);

  auto weights = std::make_shared<Value<std::vector<double>>>(
      "Weights", DataSample.Weights);

  std::shared_ptr<ComPWA::FunctionTree::TreeNode> EvaluationTree =
      std::make_shared<ComPWA::FunctionTree::TreeNode>(
          std::make_shared<Value<double>>(),
          std::make_shared<AddAll>(ParType::DOUBLE));

  auto dataTree = std::make_shared<ComPWA::FunctionTree::TreeNode>(
      std::make_shared<Value<double>>(),
      std::make_shared<MultAll>(ParType::DOUBLE));
  dataTree->addNode(createLeaf(-1));
  auto Sum = std::make_shared<TreeNode>(
      std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)));
  dataTree->addNode(Sum);
  auto WeightedLogIntensities = std::make_shared<TreeNode>(
      std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)));
  Sum->addNode(WeightedLogIntensities);
  if (weights)
    WeightedLogIntensities->addNode(createLeaf(weights));
  auto Log = std::make_shared<TreeNode>(
      std::shared_ptr<Strategy>(new LogOf(ParType::MDOUBLE)));
  WeightedLogIntensities->addNode(Log);
  Log->addNode(DataIntensityFunctionTree);

  EvaluationTree->addNode(dataTree);

  LOG(DEBUG) << "createMinLogLHEstimatorFunctionTree(): construction of LH "
                "tree finished! Performing checks ...";

  EvaluationTree->parameter();
  LOG(DEBUG) << "createMinLogLHEstimatorFunctionTree(): finished!";

  return std::make_pair(FunctionTreeEstimator(EvaluationTree, Parameters),
                        FitParList);
}

} // namespace Estimator
} // namespace ComPWA
