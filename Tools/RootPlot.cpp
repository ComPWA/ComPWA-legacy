// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "TFile.h"
#include "Tools/RootPlot.hpp"
#include "Core/ProgressBar.hpp"
#include "Tools/Integration.hpp"

using namespace ComPWA::Tools;

RootPlot::RootPlot(std::shared_ptr<ComPWA::Kinematics> kin) : Kin(kin) {}

void RootPlot::setIntensity(std::shared_ptr<ComPWA::AmpIntensity> intens) {
  PlotComponents.clear();
  PlotComponents.push_back(intens);
  ComponentNames.clear();
  ComponentNames.push_back("Intensity");
}

void RootPlot::write(std::string treePrefix, std::string fileName,
                     std::string option) {

  TFile *tf = new TFile(TString(fileName), TString(option));

  auto varNames = Kin->variableNames();
  varNames.push_back("weight");
  varNames.push_back("eff");
  
  auto varTitles = Kin->variableTitles();
  varTitles.push_back("weight");
  varTitles.push_back("#epsilon");

  size_t dataPointSize = varNames.size();
  double dataIntegral = 0.;
  // Data
  if (DataSample.size()) {
    TTree *dataTree = new TTree(TString(treePrefix + "_data"), "dataSample");
    auto t_dataSample = std::vector<double>(dataPointSize, 0.0);
    for (int i = 0; i < varNames.size(); i++)
      dataTree->Branch(TString(varNames.at(i)), &t_dataSample.at(i),
                       TString(varNames.at(i) + "/D"));

    for (auto point : DataSample) {
      // Fill branch references with dataPoint
      for (int i = 0; i < t_dataSample.size(); i++) {
        if (i < point.size())
          t_dataSample.at(i) = point.value(i);
        else if (i == point.size())
          t_dataSample.at(i) = point.weight();
        else if (i == point.size() + 1)
          t_dataSample.at(i) = point.efficiency();
        else { // Hopefully we don't arrive here
          throw std::runtime_error(
              "RootPlot::Write() | This should not happen!");
        }
      }

      dataIntegral += point.weight();
      if (point.efficiency() == 0.0) {
        LOG(ERROR) << "RootPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      dataTree->Fill();
    }
    dataTree->Write();
  }
  assert(ComponentNames.size() == PlotComponents.size());
  
  // Phase space sample
  if (PhspSample.size()) {

    double phspIntegral = std::accumulate(
        PhspSample.begin(), PhspSample.end(), 0.0,
        [&](double w, DataPoint p) { return w += p.weight(); });

    TTree *phspTree = new TTree(TString(treePrefix + "_phsp"), "phspSample");

    auto t_phspSample = std::vector<double>(dataPointSize, 0.0);
    for (int i = 0; i < varNames.size(); i++)
      phspTree->Branch(TString(varNames.at(i)), &t_phspSample.at(i),
                       TString(varNames.at(i) + "/D"));

    std::vector<double> t_weights =
        std::vector<double>(ComponentNames.size(), 0.0);
    for (int i = 0; i < ComponentNames.size(); i++)
      phspTree->Branch(TString(ComponentNames.at(i)), &t_weights.at(i),
                       TString(ComponentNames.at(i) + "/D"));

    ComPWA::ProgressBar bar(PhspSample.size());
    for (auto point : PhspSample) {
      bar.next();
      // Fill branch references with dataPoint
      for (int i = 0; i < t_phspSample.size(); i++) {
        if (i < point.size()) // variables
          t_phspSample.at(i) = point.value(i);
        else if (i == point.size()) // weight
          t_phspSample.at(i) = point.weight() * dataIntegral / phspIntegral;
        else if (i == point.size() + 1) // efficiency
          t_phspSample.at(i) = point.efficiency();
        else { // Hopefully we don't arrive here
          throw std::runtime_error(
              "RootPlot::Write() | This should not happen!");
        }
      }

      if (point.efficiency() == 0.0) {
        LOG(ERROR) << "RootPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }

      // Loop over all components that we want to plot
      for (int t = 0; t < PlotComponents.size(); t++) {
        t_weights.at(t) = PlotComponents.at(t)->intensity(point);
      }
      phspTree->Fill();
    }
    phspTree->Write();
  }
  tf->Close();
}
