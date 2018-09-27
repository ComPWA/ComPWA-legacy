// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <numeric>
#include <stdio.h>

#include "Core/DataPoint.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Core/Properties.hpp"
#include "Physics/PartialAmplitude.hpp"

#include "Tools/Plotting/ROOT/RootPlotData.hpp"

#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

namespace ComPWA {
namespace Tools {
namespace Plotting {

RootPlotData::RootPlotData(std::shared_ptr<ComPWA::Kinematics> kin,
                           std::shared_ptr<ComPWA::AmpIntensity> intens)
    : Kinematics(kin), Intensity(intens), CorrectForEfficiency(false) {}

RootPlotData::~RootPlotData() {}

/// Add sub component of the Intensity. For each event in the phase space
/// sample each component is evaluated and its value is added to the TTree.
void RootPlotData::addComponent(std::string componentName,
                                std::string intensityName, std::string title) {
  std::string ComponentLabel(title);
  if (ComponentLabel == "")
    ComponentLabel = componentName;

  std::shared_ptr<ComPWA::AmpIntensity> comp;
  try {
    comp = Intensity->component(intensityName)->component(componentName);
  } catch (std::exception &ex) {
    LOG(ERROR) << "RootPlotData::addComponent() | Component " << componentName
               << " of " << componentName << " not found in AmpIntensity "
               << Intensity->name() << ".";
    return;
  }
  AmplitudeComponents[ComponentLabel] = comp;
}

void RootPlotData::write(std::string treePrefix, std::string fileName,
                         std::string option) {

  TFile tf(fileName.c_str(), option.c_str());

  auto varNames = Kinematics->variableNames();

  // write out a final state id to name mapping
  const auto &KinProps(Kinematics->getKinematicsProperties());

  tf.mkdir("final_state_id_to_name_mapping")->cd();
  for (unsigned int i = 0; i < KinProps.FinalState.size(); ++i) {
    TParameter<int> TempPar(
        FindParticle(KinProps.ParticleList, KinProps.FinalState[i])
            .name()
            .c_str(),
        KinProps.FinalStateEventPositionMapping[i]);
    TempPar.Write();
  }
  tf.cd();

  double EventWeight(1.0);
  double EventEfficiency(1.0);

  double DataIntegral(0.0);
  //===== write data
  if (Data && Data->numEvents() > 0) {
    TTree *tree = new TTree((treePrefix + "_data").c_str(), "DataSample");
    auto DataPointValues = std::vector<double>(varNames.size(), 0.0);
    for (unsigned int i = 0; i < varNames.size(); ++i)
      tree->Branch(varNames.at(i).c_str(), &DataPointValues.at(i),
                   (varNames.at(i) + "/D").c_str());
    tree->Branch("event_weight", &EventWeight, "event_weight/D");
    tree->Branch("event_efficiency", &EventEfficiency, "event_efficiency/D");

    ComPWA::ProgressBar bar(Data->numEvents());
    for (unsigned int i = 0; i < Data->numEvents(); ++i) { // loop over data
      Event event(Data->event(i));

      EventEfficiency = 1.0;
      if (CorrectForEfficiency)
        EventEfficiency = event.efficiency();
      if (EventEfficiency == 0.0) {
        LOG(ERROR) << "DalitzPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      EventWeight = event.weight();

      /* Construct DataPoint from Event to check if dataPoint is within the
       * phase space boundaries */
      DataPoint point;
      try {
        Kinematics->convert(event, point);
      } catch (BeyondPhsp &ex) { // event outside phase, remove
        continue;
      }

      // Fill branch references with dataPoint
      for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
        DataPointValues[j] = point.value(j);
      }

      DataIntegral += point.weight();
      tree->Fill();
      bar.next();
    }
    tree->Write();
    delete (tree);
  }

  //===== write amplitude weights
  if (WeightedPhspMC && WeightedPhspMC->numEvents() > 0) {
    // calculated the total intensity integral over the phase space
    LOG(INFO) << "RootPlotData::write | calculating total intensity integral"
                 " using phase space sample...";
    double PhspIntensityIntegral(0.0);
    for (auto const &event : WeightedPhspMC->events()) { // loop over data
      DataPoint point;
      try {
        Kinematics->convert(event, point);
      } catch (BeyondPhsp &ex) { // event outside phase, remove
        continue;
      }
      PhspIntensityIntegral += event.weight() * Intensity->intensity(point);
    }

    double ScalingFactor(DataIntegral / PhspIntensityIntegral);

    TTree *tree = new TTree((treePrefix + "_weighted_phsp_MC").c_str(),
                            "WeightedPhspMCSample");
    auto DataPointValues = std::vector<double>(varNames.size(), 0.0);
    for (unsigned int i = 0; i < varNames.size(); ++i)
      tree->Branch(varNames.at(i).c_str(), &DataPointValues.at(i),
                   (varNames.at(i) + "/D").c_str());
    tree->Branch("event_weight", &EventWeight, "event_weight/D");
    tree->Branch("event_efficiency", &EventEfficiency, "event_efficiency/D");

    double IntensityWeight(0.0);
    tree->Branch("intensity", &IntensityWeight, "intensity/D");
    std::vector<double> AmplitudeComponentWeights =
        std::vector<double>(AmplitudeComponents.size(), 0.0);
    unsigned int counter(0);
    for (auto const &amp : AmplitudeComponents) {
      tree->Branch(amp.first.c_str(), &AmplitudeComponentWeights.at(counter),
                   (amp.first + "/D").c_str());
      ++counter;
    }

    ComPWA::ProgressBar bar(WeightedPhspMC->numEvents());
    for (unsigned int i = 0; i < WeightedPhspMC->numEvents();
         i++) { // loop over data
      Event event(WeightedPhspMC->event(i));

      EventEfficiency = 1.0;
      if (CorrectForEfficiency)
        EventEfficiency = event.efficiency();
      if (EventEfficiency == 0.0) {
        LOG(ERROR) << "DalitzPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      EventWeight = event.weight() * ScalingFactor;
      /* Construct DataPoint from Event to check if dataPoint is within the
       * phase space boundaries */
      DataPoint point;
      try {
        Kinematics->convert(event, point);
      } catch (BeyondPhsp &ex) { // event outside phase, remove
        continue;
      }

      // Fill branch references with dataPoint
      for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
        DataPointValues[j] = point.value(j);
      }

      IntensityWeight = Intensity->intensity(point);
      // Loop over all components that we want to plot
      counter = 0;
      for (auto const &amp : AmplitudeComponents) {
        AmplitudeComponentWeights[counter] = amp.second->intensity(point);
        ++counter;
      }

      tree->Fill();
      bar.next();
    }
    tree->Write();
    delete (tree);
  }

  //===== write data
  if (HitAndMissMC && HitAndMissMC->numEvents() > 0) {
    TTree *tree =
        new TTree((treePrefix + "_hitmiss_MC").c_str(), "Hit&MissMCSample");
    auto DataPointValues = std::vector<double>(varNames.size(), 0.0);
    for (unsigned int i = 0; i < varNames.size(); ++i)
      tree->Branch(varNames.at(i).c_str(), &DataPointValues.at(i),
                   (varNames.at(i) + "/D").c_str());
    tree->Branch("event_weight", &EventWeight, "event_weight/D");
    tree->Branch("event_efficiency", &EventEfficiency, "event_efficiency/D");

    ComPWA::ProgressBar bar(HitAndMissMC->numEvents());
    for (unsigned int i = 0; i < HitAndMissMC->numEvents();
         ++i) { // loop over data
      Event event(HitAndMissMC->event(i));

      EventEfficiency = 1.0;
      if (CorrectForEfficiency)
        EventEfficiency = event.efficiency();
      if (EventEfficiency == 0.0) {
        LOG(ERROR) << "DalitzPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      EventWeight = event.weight();

      /* Construct DataPoint from Event to check if dataPoint is within the
       * phase space boundaries */
      DataPoint point;
      try {
        Kinematics->convert(event, point);
      } catch (BeyondPhsp &ex) { // event outside phase, remove
        continue;
      }

      // Fill branch references with dataPoint
      for (unsigned int j = 0; j < DataPointValues.size(); ++j) {
        DataPointValues[j] = point.value(j);
      }

      tree->Fill();
      bar.next();
    }
    tree->Write();
    delete (tree);
  }

  tf.Close();
}
} // namespace Plotting
} // namespace Tools
} // namespace ComPWA
