// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Tools/Plotting/ROOT/RootPlotData.hpp"
#include "Core/Intensity.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Core/Properties.hpp"
#include "Data/DataTransformation.hpp"
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

void RootPlotData::writeData(
    const std::vector<DataPoint> &DataSample,
    const std::vector<std::string> &KinematicVariableNames) {
  RootFile.cd();

  double EventWeight(1.0);
  double EventEfficiency(1.0);

  double DataIntegral(0.0);

  TTree tree("data", "DataSample");
  auto DataPointValues =
      std::vector<double>(KinematicVariableNames.size(), 0.0);
  for (unsigned int i = 0; i < KinematicVariableNames.size(); ++i)
    tree.Branch(KinematicVariableNames.at(i).c_str(), &DataPointValues.at(i),
                (KinematicVariableNames.at(i) + "/D").c_str());
  tree.Branch("weight", &EventWeight, "event_weight/D");
  tree.Branch("efficiency", &EventEfficiency, "event_efficiency/D");

  ComPWA::ProgressBar bar(DataSample.size());
  for (auto const &point : DataSample) {
    EventEfficiency = point.Efficiency;
    if (EventEfficiency == 0.0) {
      LOG(ERROR) << "DalitzPlot::Fill() | Loop over "
                    "data sample: An event with zero efficiency was found! "
                    "This should not happen! We skip it!";
      continue;
    }
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
    const std::vector<DataPoint> &PhspSample,
    const std::vector<std::string> &KinematicVariableNames,
    std::shared_ptr<ComPWA::Intensity> Intensity,
    std::map<std::string, std::shared_ptr<const ComPWA::Intensity>>
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
  auto DataPointValues =
      std::vector<double>(KinematicVariableNames.size(), 0.0);
  double EventWeight(1.0);
  double EventEfficiency(1.0);
  for (unsigned int i = 0; i < KinematicVariableNames.size(); ++i)
    tree.Branch(KinematicVariableNames.at(i).c_str(), &DataPointValues.at(i),
                (KinematicVariableNames.at(i) + "/D").c_str());
  tree.Branch("weight", &EventWeight, "event_weight/D");
  tree.Branch("efficiency", &EventEfficiency, "event_efficiency/D");

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

  ComPWA::ProgressBar bar(PhspSample.size());
  for (auto const &point : PhspSample) {
    EventEfficiency = point.Efficiency;
    if (EventEfficiency == 0.0) {
      LOG(ERROR) << "DalitzPlot::Fill() | Loop over "
                    "data sample: An event with zero efficiency was found! "
                    "This should not happen! We skip it!";
      continue;
    }
    EventWeight = point.Weight;

    for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
      DataPointValues[j] = point.KinematicVariableList[j];
    }

    IntensityWeight = Intensity->evaluate(point);
    // Loop over all components that we want to plot
    counter = 0;
    for (auto const &amp : IntensityComponents) {
      AmplitudeComponentWeights[counter] = amp.second->evaluate(point);
      ++counter;
    }

    tree.Fill();
    bar.next();
  }
  tree.Write();
}

void RootPlotData::writeHitMissSample(
    const std::vector<DataPoint> &HitMissSample,
    const std::vector<std::string> &KinematicVariableNames) {

  RootFile.cd();

  TTree tree("hitmiss_data", "Hit&MissMCSample");
  auto DataPointValues =
      std::vector<double>(KinematicVariableNames.size(), 0.0);
  double EventWeight(1.0);
  double EventEfficiency(1.0);
  for (unsigned int i = 0; i < KinematicVariableNames.size(); ++i)
    tree.Branch(KinematicVariableNames.at(i).c_str(), &DataPointValues.at(i),
                (KinematicVariableNames.at(i) + "/D").c_str());
  tree.Branch("weight", &EventWeight, "event_weight/D");
  tree.Branch("efficiency", &EventEfficiency, "event_efficiency/D");

  ComPWA::ProgressBar bar(HitMissSample.size());
  for (auto const &point : HitMissSample) {
    EventEfficiency = point.Efficiency;
    if (EventEfficiency == 0.0) {
      LOG(ERROR) << "DalitzPlot::Fill() | Loop over "
                    "data sample: An event with zero efficiency was found! "
                    "This should not happen! We skip it!";
      continue;
    }
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
