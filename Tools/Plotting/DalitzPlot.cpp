// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <numeric>
#include <stdio.h>

#include "DalitzPlot.hpp"
#include "HistTools.hpp"

#include "Core/Event.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Data/DataSet.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

#include "Math/ProbFuncMathCore.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

namespace ComPWA {
namespace Tools {
namespace Plotting {

using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

void phspContour(unsigned int xsys, unsigned int ysys, unsigned int n,
                 double *xcoord, double *ycoord) {

  unsigned int num = n;
  if (num % 2 != 0) {
    num -= 1;
    LOG(INFO) << "DalitzKinematics::phspContour() | "
                 "Setting size to a even number. Assure that the size of "
                 "your arrays is "
              << num * 2 + 1 << "!";
  }
  return;
}

DalitzPlot::DalitzPlot(HelicityKinematics &kin, const std::string &name,
                       int bins)
    : Name(name), HelKin(kin), _bins(bins), _globalScale(1.0) {
  kin.createAllSubsystems();

  gStyle->SetOptStat(10); // entries only
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  gStyle->SetMarkerStyle(20);
}

void DalitzPlot::fill(const std::vector<ComPWA::Event> &data, bool normalize,
                      const std::string &name, const std::string &title,
                      Color_t color) {

  Data::DataSet dataset = HelKin.convert(data);

  DalitzHisto hist(HelKin, name, title, _bins, color);
  hist.setStats(0);

  hist.fill(dataset);

  if (normalize)
    _globalScale = hist.integral();

  _plotHistograms.push_back(std::move(hist));
}

void DalitzPlot::fill(const std::vector<ComPWA::Event> &data, Intensity &intens,
                      bool normalize, const std::string &name,
                      const std::string &title, Color_t color) {
  Data::DataSet dataset = HelKin.convert(data);
  DalitzHisto hist(HelKin, name, title, _bins, color);
  hist.setStats(0);
  auto Intensities = intens.evaluate(dataset.Data);
  hist.fill(dataset, Intensities);

  if (normalize)
    _globalScale = hist.integral();

  _plotHistograms.push_back(std::move(hist));
}

void DalitzPlot::plot() {
  for (auto &pl : _plotHistograms)
    pl.scale(_globalScale / pl.integral());

  //=== generate contour
  double xpoints[4001], ypoints[4001];
  phspContour(0, 1, 2000, xpoints, ypoints);
  TGraph m23m13_contour(4001, xpoints, ypoints);
  m23m13_contour.SetMarkerStyle(1);
  m23m13_contour.SetLineColor(kRed);
  m23m13_contour.SetMarkerColor(kRed);
  m23m13_contour.SetMarkerSize(0.0);
  m23m13_contour.SetTitle("phspContour");
  m23m13_contour.SetFillColor(kWhite);
  phspContour(0, 2, 2000, xpoints, ypoints);
  TGraph m23m12_contour(4001, xpoints, ypoints);
  m23m12_contour.SetMarkerStyle(1);
  m23m12_contour.SetLineColor(kRed);
  m23m12_contour.SetMarkerColor(kRed);
  m23m12_contour.SetMarkerSize(0.0);
  m23m12_contour.SetTitle("phspContour");
  m23m12_contour.SetFillColor(kWhite);
  phspContour(2, 1, 2000, xpoints, ypoints);
  TGraph m12m13_contour(4001, xpoints, ypoints);
  m12m13_contour.SetMarkerStyle(1);
  m12m13_contour.SetLineColor(kRed);
  m12m13_contour.SetMarkerColor(kRed);
  m12m13_contour.SetMarkerSize(0.0);
  m12m13_contour.SetTitle("phspContour");
  m12m13_contour.SetFillColor(kWhite);

  //----- plotting invariant mass distributions -----
  TCanvas *c2 = new TCanvas("invmass", "invmass", 50, 50, 2400, 800);
  c2->Divide(3, 1);
  c2->cd(1);
  CreateHist("mSq_(1,2)"); // Plotting mKKsq
  c2->cd(2);
  CreateHist("mSq_(0,2)"); // Plotting mKSK+sq
  c2->cd(3);
  CreateHist("mSq_(0,1)"); // Plotting mKSK+sq

  //----- Write to TFile -----
  TFile *tf2 = new TFile(Name + ".root", "recreate");
  if (tf2->IsZombie()) {
    std::cout << "Error opening output file" << std::endl;
    exit(-1);
  }
  m12m13_contour.Write("m12m13_contour", TObject::kOverwrite, 0);
  m23m12_contour.Write("m23m12_contour", TObject::kOverwrite, 0);
  m23m13_contour.Write("m23m13_contour", TObject::kOverwrite, 0);
  c2->Write("invmass", TObject::kOverwrite, 0);

  // Save data trees and histograms
  tf2->mkdir("hist");
  tf2->cd("hist");
  for (auto &pl : _plotHistograms)
    pl.write();

  // Write some canvas to single files
  c2->Print(Name + "-invmass.pdf");

  tf2->Close();

  return;
}

void DalitzPlot::CreateHist(std::string Name) {
  if (!_plotHistograms.size())
    return;

  std::vector<TH1D *> v;
  std::vector<TString> options;
  v.push_back(_plotHistograms.at(0).getHistogram(Name));
  options.push_back("E1");
  for (unsigned int t = 1; t < _plotHistograms.size(); ++t) {
    v.push_back(_plotHistograms.at(t).getHistogram(Name));
    options.push_back("Sames,Hist");
  }

  drawPull(v, options);
}

//===================== DalitzHisto =====================
DalitzHisto::DalitzHisto(HelicityKinematics &helkin, std::string name,
                         std::string title, unsigned int bins, Color_t color)
    : Name(name), Title(title), NumBins(bins), Integral(0.0) {

  // Initialize TTree
  Tree = std::unique_ptr<TTree>(new TTree(TString(Name), TString(Title)));

  // Adding branches to TTree
  Tree->Branch("sample", &BranchPoint);
  Tree->Branch("efficiency", &BranchEff, "eff/D");
  Tree->Branch("weight", &BranchWeight, "weight/D");

  char label[60];

  // mass23sq
  auto sys23(helkin.registerSubSystem({1}, {2}, {0}, {}));
  auto m23sq_limit = helkin.getInvariantMassBounds(std::get<0>(sys23));
  double m23sq_min = m23sq_limit.first;
  double m23sq_max = m23sq_limit.second;
  TH1D hist = TH1D(TString(Name + "m23sq"), TString(Title), NumBins, m23sq_min,
                   m23sq_max);
  double binWidth = (double)(m23sq_min - m23sq_max) / NumBins;
  sprintf(label, "Entries /%f.3", binWidth);
  hist.GetYaxis()->SetTitle("# [" + TString(label) + "]");
  hist.GetXaxis()->SetTitle("m_{23}^{2} [GeV/c^{2}]");
  hist.Sumw2();
  Hists1D.insert(std::make_pair(std::get<0>(sys23), hist));

  // mass13sq
  auto sys13(helkin.registerSubSystem({0}, {2}, {1}, {}));
  auto m13sq_limit = helkin.getInvariantMassBounds(std::get<0>(sys13));
  double m13sq_min = m13sq_limit.first;
  double m13sq_max = m13sq_limit.second;
  hist = TH1D(TString(Name + "m13sq"), TString(Title), NumBins, m13sq_min,
              m13sq_max);
  binWidth = (double)(m13sq_min - m13sq_max) / NumBins;
  sprintf(label, "Entries /%f.3", binWidth);
  hist.GetYaxis()->SetTitle("# [" + TString(label) + "]");
  hist.GetXaxis()->SetTitle("m_{13}^{2} [GeV/c^{2}]");
  hist.Sumw2();
  Hists1D.insert(std::make_pair(std::get<0>(sys13), hist));

  // mass12sq
  auto sys12(helkin.registerSubSystem({0}, {1}, {2}, {}));
  auto m12sq_limit = helkin.getInvariantMassBounds(std::get<0>(sys12));
  double m12sq_min = m12sq_limit.first;
  double m12sq_max = m12sq_limit.second;
  hist = TH1D(TString(Name + "m12sq"), TString(Title), NumBins, m12sq_min,
              m12sq_max);
  binWidth = (double)(m12sq_min - m12sq_max) / NumBins;
  sprintf(label, "Entries /%f.3", binWidth);
  hist.GetYaxis()->SetTitle("# [" + TString(label) + "]");
  hist.GetXaxis()->SetTitle("m_{12}^{2} [GeV/c^{2}]");
  hist.Sumw2();
  Hists1D.insert(std::make_pair(std::get<0>(sys12), hist));

  TH2D hist2d = TH2D(TString(name + "_m23sqm13sq"), TString(title), NumBins,
                     m23sq_min, m23sq_max, NumBins, m13sq_min, m13sq_max);
  hist2d.GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
  hist2d.GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
  Hists2D.insert(std::make_pair(
      std::make_pair(std::get<0>(sys23), std::get<0>(sys13)), hist2d));

  hist2d = TH2D(TString(name + "_m23sqm12sq"), TString(title), NumBins,
                m23sq_min, m23sq_max, NumBins, m12sq_min, m12sq_max);
  hist2d.GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}");
  hist2d.GetYaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
  Hists2D.insert(std::make_pair(
      std::make_pair(std::get<0>(sys23), std::get<0>(sys12)), hist2d));

  hist2d = TH2D(TString(name + "_m12sqm13sq"), TString(title), NumBins,
                m12sq_min, m12sq_max, NumBins, m13sq_min, m13sq_max);
  hist2d.GetXaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
  hist2d.GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
  Hists2D.insert(std::make_pair(
      std::make_pair(std::get<0>(sys12), std::get<0>(sys13)), hist2d));

  hist2d = TH2D(TString(name + "_m23sqCosTheta"), TString(title), NumBins,
                m23sq_min, m23sq_max, NumBins, -1, 1);
  hist2d.GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
  hist2d.GetYaxis()->SetTitle("#cos(#Theta)_{KK}");
  Hists2D.insert(std::make_pair(
      std::make_pair(std::get<0>(sys12), std::get<0>(sys12)), hist2d));

  for (auto x : Hists2D) {
    x.second.GetXaxis()->SetNdivisions(508);
    x.second.GetZaxis()->SetTitle("Entries");
  }

  setColor(color);
  return;
}

void DalitzHisto::fill(const Data::DataSet &sample,
                       std::vector<double> Intensities) {
  if (Intensities.size() == 0) {
    Intensities = std::vector<double>(sample.Weights.size(), 1.0);
  }
  if (!sample.Data.size())
    throw std::runtime_error("DalitzHist::fill() | Empty data sample.");
  if (Intensities.size() != sample.Weights.size())
    throw std::runtime_error("DalitzHist::fill() | Vector of weights and "
                             "sample do not have the same length");

  double weight =
      std::accumulate(sample.Weights.begin(), sample.Weights.end(), 0.0);

  Integral += weight;

  std::vector<double> weights;
  for (size_t i = 0; i < Intensities.size(); ++i) {
    weights.push_back(sample.Weights[i] * Intensities[i]);
  }

  for (auto h1 : Hists1D) {
    auto const &data = sample.Data.at(h1.first);
    for (size_t i = 0; i < data.size(); ++i) {
      h1.second.Fill(data[i], weights[i]);
    }
  }

  for (auto h2 : Hists2D) {
    auto const &data1 = sample.Data.at(h2.first.first);
    auto const &data2 = sample.Data.at(h2.first.second);
    for (size_t i = 0; i < data1.size(); ++i) {
      h2.second.Fill(data1[i], data2[i], weights[i]);
    }
  }

  auto const &m23sq = sample.Data.at("mSq_(1,2)");
  auto const &m13sq = sample.Data.at("mSq_(0,2)");
  auto const &m12sq = sample.Data.at("mSq_(0,1)");
  for (size_t i = 0; i < weights.size(); ++i) {
    BranchPoint = {m23sq[i], m13sq[i], m12sq[i]};
    BranchWeight = weights[i];
    BranchEff = 1.0;
    Tree->Fill();
  }
}

void DalitzHisto::setStats(bool b) {
  for (auto x : Hists1D) {
    x.second.SetStats(b);
  }
  for (auto x : Hists2D) {
    x.second.SetStats(b);
  }
}

void DalitzHisto::scale(double w) {
  for (auto x : Hists1D) {
    x.second.Scale(w);
  }
  for (auto x : Hists2D) {
    x.second.Scale(w);
  }
}

void DalitzHisto::setColor(Color_t color) {
  for (auto x : Hists1D) {
    x.second.SetLineColor(color);
    x.second.SetMarkerColor(color);
  }
}

TH1D *DalitzHisto::getHistogram(std::string Name) { return &Hists1D.at(Name); }

TH2D *DalitzHisto::getHistogram2D(std::pair<std::string, std::string> Names) {
  return &Hists2D.at(Names);
}

void DalitzHisto::write() {
  Tree->Write(TString(Name) + "Tree");
  gDirectory->mkdir(TString(Name) + "_hist");
  gDirectory->cd(TString(Name) + "_hist");
  for (auto const &x : Hists1D) {
    x.second.Write();
  }
  for (auto const &x : Hists2D) {
    x.second.Write();
  }
  gDirectory->cd("../");
}

} // namespace Plotting
} // namespace Tools
} // namespace ComPWA
