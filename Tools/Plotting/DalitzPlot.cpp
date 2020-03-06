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

void phspContour(unsigned int SysX, unsigned int SysY, unsigned int n,
                 double *CoordX, double *CoordY) {

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

DalitzPlot::DalitzPlot(HelicityKinematics &Kinematics, const std::string &Name,
                       int Bins)
    : Name(Name), Kinematics(Kinematics), Bins(Bins), GlobalScale(1.0) {
  Kinematics.createAllSubsystems();

  gStyle->SetOptStat(10); // entries only
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  gStyle->SetMarkerStyle(20);
}

void DalitzPlot::fill(const ComPWA::EventCollection &Data, bool Normalize,
                      const std::string &Name, const std::string &Title,
                      Color_t Color) {

  Data::DataSet DataSet = Kinematics.convert(Data);

  DalitzHisto Histogram(Kinematics, Name, Title, Bins, Color);
  Histogram.setStats(0);

  Histogram.fill(DataSet);

  if (Normalize)
    GlobalScale = Histogram.integral();

  PlotHistograms.push_back(std::move(Histogram));
}

void DalitzPlot::fill(const ComPWA::EventCollection &Data, Intensity &Intensity,
                      bool Normalize, const std::string &Name,
                      const std::string &Title, Color_t Color) {
  Data::DataSet DataSet = Kinematics.convert(Data);
  DalitzHisto Histogram(Kinematics, Name, Title, Bins, Color);
  Histogram.setStats(0);
  auto Intensities = Intensity.evaluate(DataSet.Data);
  Histogram.fill(DataSet, Intensities);

  if (Normalize)
    GlobalScale = Histogram.integral();

  PlotHistograms.push_back(std::move(Histogram));
}

void DalitzPlot::plot() {
  for (auto &Plot : PlotHistograms)
    Plot.scale(GlobalScale / Plot.integral());

  //=== generate contour
  double PointsX[4001], PointsY[4001];
  phspContour(0, 1, 2000, PointsX, PointsY);
  TGraph ContourM23M13(4001, PointsX, PointsY);
  ContourM23M13.SetMarkerStyle(1);
  ContourM23M13.SetLineColor(kRed);
  ContourM23M13.SetMarkerColor(kRed);
  ContourM23M13.SetMarkerSize(0.0);
  ContourM23M13.SetTitle("phspContour");
  ContourM23M13.SetFillColor(kWhite);
  phspContour(0, 2, 2000, PointsX, PointsY);
  TGraph ContourM23M12(4001, PointsX, PointsY);
  ContourM23M12.SetMarkerStyle(1);
  ContourM23M12.SetLineColor(kRed);
  ContourM23M12.SetMarkerColor(kRed);
  ContourM23M12.SetMarkerSize(0.0);
  ContourM23M12.SetTitle("phspContour");
  ContourM23M12.SetFillColor(kWhite);
  phspContour(2, 1, 2000, PointsX, PointsY);
  TGraph ContourM12M13(4001, PointsX, PointsY);
  ContourM12M13.SetMarkerStyle(1);
  ContourM12M13.SetLineColor(kRed);
  ContourM12M13.SetMarkerColor(kRed);
  ContourM12M13.SetMarkerSize(0.0);
  ContourM12M13.SetTitle("phspContour");
  ContourM12M13.SetFillColor(kWhite);

  //----- plotting invariant mass distributions -----
  TCanvas *Canvas = new TCanvas("invmass", "invmass", 50, 50, 2400, 800);
  Canvas->Divide(3, 1);
  Canvas->cd(1);
  CreateHist("mSq_(1,2)"); // Plotting mKKsq
  Canvas->cd(2);
  CreateHist("mSq_(0,2)"); // Plotting mKSK+sq
  Canvas->cd(3);
  CreateHist("mSq_(0,1)"); // Plotting mKSK+sq

  //----- Write to TFile -----
  TFile *File = new TFile(Name + ".root", "recreate");
  if (File->IsZombie()) {
    std::cout << "Error opening output file" << std::endl;
    exit(-1);
  }
  ContourM12M13.Write("m12m13_contour", TObject::kOverwrite, 0);
  ContourM23M12.Write("m23m12_contour", TObject::kOverwrite, 0);
  ContourM23M13.Write("m23m13_contour", TObject::kOverwrite, 0);
  Canvas->Write("invmass", TObject::kOverwrite, 0);

  // Save data trees and histograms
  File->mkdir("hist");
  File->cd("hist");
  for (auto &Plot : PlotHistograms)
    Plot.write();

  // Write some canvas to single files
  Canvas->Print(Name + "-invmass.pdf");

  File->Close();
}

void DalitzPlot::CreateHist(std::string Name) {
  if (!PlotHistograms.size())
    return;

  std::vector<TH1D *> Histograms;
  std::vector<TString> Options;
  Histograms.push_back(PlotHistograms.at(0).getHistogram(Name));
  Options.push_back("E1");
  for (unsigned int t = 1; t < PlotHistograms.size(); ++t) {
    Histograms.push_back(PlotHistograms.at(t).getHistogram(Name));
    Options.push_back("Sames,Hist");
  }

  drawPull(Histograms, Options);
}

//===================== DalitzHisto =====================
DalitzHisto::DalitzHisto(HelicityKinematics &Kinematics, std::string Name,
                         std::string Title, unsigned int Bins, Color_t Color)
    : Name(Name), Title(Title), NumBins(Bins), Integral(0.0) {

  // Initialize TTree
  Tree = std::unique_ptr<TTree>(new TTree(TString(Name), TString(Title)));

  // Adding branches to TTree
  Tree->Branch("sample", &BranchPoint);
  Tree->Branch("efficiency", &BranchEff, "eff/D");
  Tree->Branch("weight", &BranchWeight, "weight/D");

  char label[60];

  // mass23sq
  auto sys23(Kinematics.registerSubSystem({1}, {2}, {0}, {}));
  auto m23sq_limit = Kinematics.getInvariantMassBounds(std::get<0>(sys23));
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
  auto sys13(Kinematics.registerSubSystem({0}, {2}, {1}, {}));
  auto m13sq_limit = Kinematics.getInvariantMassBounds(std::get<0>(sys13));
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
  auto sys12(Kinematics.registerSubSystem({0}, {1}, {2}, {}));
  auto m12sq_limit = Kinematics.getInvariantMassBounds(std::get<0>(sys12));
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

  TH2D hist2d = TH2D(TString(Name + "_m23sqm13sq"), TString(Title), NumBins,
                     m23sq_min, m23sq_max, NumBins, m13sq_min, m13sq_max);
  hist2d.GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
  hist2d.GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
  Hists2D.insert(std::make_pair(
      std::make_pair(std::get<0>(sys23), std::get<0>(sys13)), hist2d));

  hist2d = TH2D(TString(Name + "_m23sqm12sq"), TString(Title), NumBins,
                m23sq_min, m23sq_max, NumBins, m12sq_min, m12sq_max);
  hist2d.GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}");
  hist2d.GetYaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
  Hists2D.insert(std::make_pair(
      std::make_pair(std::get<0>(sys23), std::get<0>(sys12)), hist2d));

  hist2d = TH2D(TString(Name + "_m12sqm13sq"), TString(Title), NumBins,
                m12sq_min, m12sq_max, NumBins, m13sq_min, m13sq_max);
  hist2d.GetXaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
  hist2d.GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
  Hists2D.insert(std::make_pair(
      std::make_pair(std::get<0>(sys12), std::get<0>(sys13)), hist2d));

  hist2d = TH2D(TString(Name + "_m23sqCosTheta"), TString(Title), NumBins,
                m23sq_min, m23sq_max, NumBins, -1, 1);
  hist2d.GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
  hist2d.GetYaxis()->SetTitle("#cos(#Theta)_{KK}");
  Hists2D.insert(std::make_pair(
      std::make_pair(std::get<0>(sys12), std::get<0>(sys12)), hist2d));

  for (auto &x : Hists2D) {
    x.second.GetXaxis()->SetNdivisions(508);
    x.second.GetZaxis()->SetTitle("Entries");
  }

  setColor(Color);
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

  for (auto &h1 : Hists1D) {
    auto const &data = sample.Data.at(h1.first);
    for (size_t i = 0; i < data.size(); ++i) {
      h1.second.Fill(data[i], weights[i]);
    }
  }

  for (auto &h2 : Hists2D) {
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
  for (auto &x : Hists1D) {
    x.second.SetStats(b);
  }
  for (auto &x : Hists2D) {
    x.second.SetStats(b);
  }
}

void DalitzHisto::scale(double w) {
  for (auto &x : Hists1D) {
    x.second.Scale(w);
  }
  for (auto &x : Hists2D) {
    x.second.Scale(w);
  }
}

void DalitzHisto::setColor(Color_t color) {
  for (auto &x : Hists1D) {
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
