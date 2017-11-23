

// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "TFile.h"
#include "Tools/RootPlot.hpp"
#include "Core/ProgressBar.hpp"
#include "Tools/Integration.hpp"

using namespace ComPWA::Tools;

RootPlot::RootPlot(std::shared_ptr<ComPWA::Kinematics> kin) : kin_(kin) {}

void RootPlot::SetFitAmp(std::shared_ptr<ComPWA::AmpIntensity> intens) {
  _plotComponents.clear();
  _plotComponents.push_back(intens);
  _componentNames.clear();
  _componentNames.push_back("Intensity");
}

void RootPlot::Write(std::string treePrefix, std::string fileName,
                     std::string option) {

  TFile *tf = new TFile(TString(fileName), TString(option));

  auto varNames = kin_->GetVarNames();
  varNames.push_back("weight");
  varNames.push_back("eff");
  
  auto varTitles = kin_->GetVarTitles();
  varTitles.push_back("weight");
  varTitles.push_back("#epsilon");

  size_t dataPointSize = varNames.size();
  double dataIntegral = 0.;
  // Data
  if (s_data.size()) {
    TTree *dataTree = new TTree(TString(treePrefix + "_data"), "dataSample");
    auto t_dataSample = std::vector<double>(dataPointSize, 0.0);
    for (int i = 0; i < varNames.size(); i++)
      dataTree->Branch(TString(varNames.at(i)), &t_dataSample.at(i),
                       TString(varNames.at(i) + "/D"));

    for (auto point : s_data) {
      // Fill branch references with dataPoint
      for (int i = 0; i < t_dataSample.size(); i++) {
        if (i < point.Size())
          t_dataSample.at(i) = point.GetValue(i);
        else if (i == point.Size())
          t_dataSample.at(i) = point.GetWeight();
        else if (i == point.Size() + 1)
          t_dataSample.at(i) = point.GetEfficiency();
        else { // Hopefully we don't arrive here
          throw std::runtime_error(
              "RootPlot::Write() | This should not happen!");
        }
      }

      dataIntegral += point.GetWeight();
      if (point.GetEfficiency() == 0.0) {
        LOG(error) << "RootPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      dataTree->Fill();
    }
    dataTree->Write();
  }

  // Phase space sample
  if (s_phsp.size()) {

    double phspIntegral = std::accumulate(
        s_phsp.begin(), s_phsp.end(), 0.0,
        [&](double w, dataPoint p) { return w += p.GetWeight(); });

    TTree *phspTree = new TTree(TString(treePrefix + "_phsp"),
                                "phspSample including amplitude weights");

    auto t_phspSample = std::vector<double>(dataPointSize, 0.0);
    for (int i = 0; i < varNames.size(); i++)
      phspTree->Branch(TString(varNames.at(i)), &t_phspSample.at(i),
                       TString(varNames.at(i) + "/D"));

    std::vector<double> t_weights =
        std::vector<double>(_componentNames.size(), 0.0);
    for (int i = 0; i < _componentNames.size(); i++)
      phspTree->Branch(TString(_componentNames.at(i)), &t_weights.at(i),
                       TString(_componentNames.at(i) + "/D"));

    ComPWA::progressBar bar(s_phsp.size());
    for (auto point : s_phsp) {
      bar.nextEvent();
      // Fill branch references with dataPoint
      for (int i = 0; i < t_phspSample.size(); i++) {
        if (i < point.Size())
          t_phspSample.at(i) = point.GetValue(i);
        else if (i == point.Size())
          t_phspSample.at(i) = point.GetWeight() * dataIntegral / phspIntegral;
        else if (i == point.Size() + 1)
          t_phspSample.at(i) = point.GetEfficiency();
        else { // Hopefully we don't arrive here
          throw std::runtime_error(
              "RootPlot::Write() | This should not happen!");
        }
      }

      if (point.GetEfficiency() == 0.0) {
        LOG(error) << "RootPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }

      // Loop over all components that we want to plot
      for (int t = 0; t < _plotComponents.size(); t++) {
        t_weights.at(t) = _plotComponents.at(t)->Intensity(point);
      }
      phspTree->Fill();
    }
    phspTree->Write();
  }
  tf->Close();
}
