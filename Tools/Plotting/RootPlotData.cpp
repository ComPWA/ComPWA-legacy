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

void RootPlotData::createDirectory(std::string Name) {
  RootFile.cd();
  RootFile.mkdir(Name.c_str());
  RootFile.cd(Name.c_str());
}

void RootPlotData::writeData(const Data::DataSet &DataSample,
                             std::string TreeName) {
  double EventWeight(1.0);

  double DataIntegral(0.0);

  TTree tree(TreeName.c_str(), "DataSample");
  auto const &KinematicVariableNames = DataSample.VariableNames;
  auto DataPointValues =
      std::vector<double>(KinematicVariableNames.size(), 0.0);
  for (unsigned int i = 0; i < KinematicVariableNames.size(); ++i)
    tree.Branch(KinematicVariableNames.at(i).c_str(), &DataPointValues.at(i),
                (KinematicVariableNames.at(i) + "/D").c_str());
  tree.Branch("weight", &EventWeight, "event_weight/D");

  ComPWA::ProgressBar bar(DataSample.Weights.size());
  for (size_t i = 0; i < DataSample.Weights.size(); ++i) {
    EventWeight = DataSample.Weights[i];

    for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
      DataPointValues[j] = DataSample.Data[j][i];
    }

    DataIntegral += DataSample.Weights[i];
    tree.Fill();
    bar.next();
  }
  tree.Write();
}

void RootPlotData::writeIntensityWeightedPhspSample(
    const Data::DataSet &PhspSample, ComPWA::Intensity &Intensity,
    std::string TreeName,
    std::map<std::string, std::shared_ptr<ComPWA::Intensity>>
        IntensityComponents) {

  LOG(INFO) << "RootPlotData::write | calculating total intensity integral"
               " using phase space sample...";

  // double PhspIntensityIntegral(0.0);
  // for (auto const &dp : PhspSample.DataPoints) { // loop over data
  //  PhspIntensityIntegral += dp.Weight * Intensity->evaluate(dp);
  //}

  //   double ScalingFactor(DataIntegral / PhspIntensityIntegral);

  TTree tree(TreeName.c_str(), "WeightedPhspMCSample");
  auto const &KinematicVariableNames = PhspSample.VariableNames;
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

  auto Intensities = Intensity.evaluate(PhspSample.Data);
  std::vector<std::vector<double>> Components;
  for (auto amp : IntensityComponents) {
    Components.push_back(amp.second->evaluate(PhspSample.Data));
  }

  ComPWA::ProgressBar bar(PhspSample.Weights.size());
  for (size_t i = 0; i < PhspSample.Weights.size(); ++i) {
    EventWeight = PhspSample.Weights[i];

    for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
      DataPointValues[j] = PhspSample.Data[j][i];
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

void RootPlotData::writeHitMissSample(const Data::DataSet &HitMissSample,
                                      std::string TreeName) {
  TTree tree(TreeName.c_str(), "Hit&MissMCSample");
  auto const &KinematicVariableNames = HitMissSample.VariableNames;
  auto DataPointValues =
      std::vector<double>(KinematicVariableNames.size(), 0.0);
  double EventWeight(1.0);

  for (unsigned int i = 0; i < KinematicVariableNames.size(); ++i)
    tree.Branch(KinematicVariableNames.at(i).c_str(), &DataPointValues.at(i),
                (KinematicVariableNames.at(i) + "/D").c_str());
  tree.Branch("weight", &EventWeight, "event_weight/D");

  ComPWA::ProgressBar bar(HitMissSample.Weights.size());
  for (size_t i = 0; i < HitMissSample.Weights.size(); ++i) {
    EventWeight = HitMissSample.Weights[i];

    for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
      DataPointValues[j] = HitMissSample.Data[j][i];
    }

    tree.Fill();
    bar.next();
  }
  tree.Write();
}

} // namespace Plotting
} // namespace Tools
} // namespace ComPWA
