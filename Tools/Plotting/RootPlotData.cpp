// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Tools/Plotting/RootPlotData.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"

#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"

namespace ComPWA {
namespace Tools {
namespace Plotting {

RootPlotData::RootPlotData(
    const Physics::ParticleStateTransitionKinematicsInfo &KinematicsInfo,
    const std::string &filename, const std::string &option)
    : RootFile(filename.c_str(), option.c_str()) {

  // write final state id to name mapping to file

  RootFile.cd();
  RootFile.mkdir("final_state_id_to_name_mapping")->cd();
  for (auto const x : KinematicsInfo.getFinalStateIDToNameMapping()) {
    TParameter<int> TempPar(x.second.c_str(), x.first);
    TempPar.Write();
  }
  RootFile.cd();
}

void RootPlotData::writeData(const Data::DataSet &DataSample) {
  RootFile.cd();

  double EventWeight(1.0);

  double DataIntegral(0.0);

  TTree tree("data", "DataSample");
  auto const &KinematicVariableNames = DataSample.getKinematicVariableNames();
  auto DataPointValues =
      std::vector<double>(KinematicVariableNames.size(), 0.0);
  for (unsigned int i = 0; i < KinematicVariableNames.size(); ++i)
    tree.Branch(KinematicVariableNames.at(i).c_str(), &DataPointValues.at(i),
                (KinematicVariableNames.at(i) + "/D").c_str());
  tree.Branch("weight", &EventWeight, "event_weight/D");

  auto const &DataPoints = DataSample.getDataPointList();
  ComPWA::ProgressBar bar(DataPoints.size());
  for (auto const &point : DataPoints) {
    EventWeight = point.Weight;

    for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
      DataPointValues[j] = point.KinematicVariableList[j];
    }

    DataIntegral += point.Weight;
    tree.Fill();
    bar.next();
  }
  tree.Write();
}

void RootPlotData::writeIntensityWeightedPhspSample(
    const Data::DataSet &PhspSample,
    std::shared_ptr<ComPWA::Intensity> Intensity,
    std::map<std::string, std::shared_ptr<ComPWA::Intensity>>
        IntensityComponents) {

  if (!Intensity) {
    throw std::runtime_error("RootPlotData::write: Supplied a phase space "
                             "sample, but no intensity for weighting!");
  }

  LOG(INFO) << "RootPlotData::write | calculating total intensity integral"
               " using phase space sample...";

  RootFile.cd();

  // double PhspIntensityIntegral(0.0);
  // for (auto const &dp : PhspSample.DataPoints) { // loop over data
  //  PhspIntensityIntegral += dp.Weight * Intensity->evaluate(dp);
  //}

  //   double ScalingFactor(DataIntegral / PhspIntensityIntegral);

  TTree tree("intensity_weighted_phspdata", "WeightedPhspMCSample");
  auto const &KinematicVariableNames = PhspSample.getKinematicVariableNames();
  auto DataPointValues =
      std::vector<double>(KinematicVariableNames.size(), 0.0);
  double EventWeight(1.0);

  for (unsigned int i = 0; i < KinematicVariableNames.size(); ++i)
    tree.Branch(KinematicVariableNames.at(i).c_str(), &DataPointValues.at(i),
                (KinematicVariableNames.at(i) + "/D").c_str());
  tree.Branch("weight", &EventWeight, "event_weight/D");

  double IntensityWeight(0.0);
  tree.Branch("intensity", &IntensityWeight, "intensity/D");
  std::vector<double> AmplitudeComponentWeights =
      std::vector<double>(IntensityComponents.size(), 0.0);
  unsigned int counter(0);
  for (auto const &amp : IntensityComponents) {
    tree.Branch(amp.first.c_str(), &AmplitudeComponentWeights.at(counter),
                (amp.first + "/D").c_str());
    ++counter;
  }

  auto const &PhspPoints = PhspSample.getParameterList();
  std::vector<std::vector<double>> PhspDataPoints;
  for (size_t i = 0; i < PhspPoints.mDoubleValues().size() - 1; ++i) {
    PhspDataPoints.push_back(PhspPoints.mDoubleValue(i)->values());
  }
  auto PhspDataPointWeights =
      PhspPoints.mDoubleValue(PhspPoints.mDoubleValues().size() - 1)->values();

  auto Intensities = Intensity->evaluate(PhspDataPoints);
  std::vector<std::vector<double>> Components;
  for (auto amp : IntensityComponents) {
    Components.push_back(amp.second->evaluate(PhspDataPoints));
  }

  ComPWA::ProgressBar bar(PhspDataPointWeights.size());
  for (size_t i = 0; i < PhspDataPointWeights.size(); ++i) {
    EventWeight = PhspDataPointWeights[i];

    for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
      DataPointValues[j] = PhspDataPoints[j][i];
    }
    IntensityWeight = Intensities[i];
    for (size_t j = 0; j < Components.size(); ++j) {
      AmplitudeComponentWeights[j] = Components[j][i];
    }

    tree.Fill();
    bar.next();
  }
  tree.Write();
}

void RootPlotData::writeHitMissSample(const Data::DataSet &HitMissSample) {

  RootFile.cd();

  TTree tree("hitmiss_data", "Hit&MissMCSample");
  auto const &KinematicVariableNames =
      HitMissSample.getKinematicVariableNames();
  auto DataPointValues =
      std::vector<double>(KinematicVariableNames.size(), 0.0);
  double EventWeight(1.0);

  for (unsigned int i = 0; i < KinematicVariableNames.size(); ++i)
    tree.Branch(KinematicVariableNames.at(i).c_str(), &DataPointValues.at(i),
                (KinematicVariableNames.at(i) + "/D").c_str());
  tree.Branch("weight", &EventWeight, "event_weight/D");

  auto const &DataPoints = HitMissSample.getDataPointList();
  ComPWA::ProgressBar bar(DataPoints.size());
  for (auto const &point : DataPoints) {
    EventWeight = point.Weight;

    for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
      DataPointValues[j] = point.KinematicVariableList[j];
    }

    tree.Fill();
    bar.next();
  }
  tree.Write();
}

} // namespace Plotting
} // namespace Tools
} // namespace ComPWA
