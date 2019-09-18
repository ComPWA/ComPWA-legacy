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

DalitzPlot::DalitzPlot(HelicityKinematics &kin, std::string name, int bins)
    : Name(name), HelKin(kin), _bins(bins), _globalScale(1.0){
      
  gStyle->SetOptStat(10); // entries only
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  gStyle->SetMarkerStyle(20);
}

void DalitzPlot::fill(const std::vector<ComPWA::Event> &data, bool normalize,
                      std::string name, std::string title, Color_t color) {

  Data::DataSet dataset = Data::convertEventsToDataSet(data, HelKin);
  
  DalitzHisto hist(HelKin, name, title, _bins, color);
  hist.setStats(0);

  hist.fill(HelKin, dataset);

  if(normalize)
    _globalScale = hist.integral();

  _plotHistograms.push_back(std::move(hist));
}

void DalitzPlot::fill(const std::vector<ComPWA::Event> &data,
                      FunctionTreeIntensity &intens, bool normalize,
                      std::string name, std::string title, Color_t color) {

  Data::DataSet dataset = Data::convertEventsToDataSet(data, HelKin);
  DalitzHisto hist(HelKin, name, title, _bins, color);
  hist.setStats(0);
  auto Intensities = intens.evaluate(dataset.Data);
  hist.fill(HelKin, dataset, Intensities);

  if(normalize)
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
  CreateHist(0); // Plotting mKKsq
  c2->cd(2);
  CreateHist(1); // Plotting mKSK+sq
  c2->cd(3);
  CreateHist(2); // Plotting mKSK+sq

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
  for(auto &pl:_plotHistograms)
    pl.write();

  // Write some canvas to single files
  c2->Print(Name + "-invmass.pdf");

  tf2->Close();

  return;
}

void DalitzPlot::CreateHist(unsigned int id) {
  if (!_plotHistograms.size())
    return;
  
  std::vector<TH1D *> v;
  std::vector<TString> options;
  v.push_back(_plotHistograms.at(0).getHistogram(id));
  options.push_back("E1");
  for (unsigned int t = 1; t < _plotHistograms.size(); ++t) {
    v.push_back(_plotHistograms.at(t).getHistogram(id));
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
  unsigned int sys23(helkin.addSubSystem({1}, {2}, {0}, {}));
  auto m23sq_limit = helkin.invMassBounds(sys23);
  double m23sq_min = m23sq_limit.first;
  double m23sq_max = m23sq_limit.second;
  Arr.push_back(
      TH1D("m23sq", "m_{23}^{2} [GeV/c^{2}]", NumBins, m23sq_min, m23sq_max));
  double binWidth = (double)(m23sq_min - m23sq_max) / NumBins;
  sprintf(label, "Entries /%f.3", binWidth);
  Arr.back().GetYaxis()->SetTitle("# [" + TString(label) + "]");
  Arr.back().GetXaxis()->SetTitle("m_{23}^{2} [GeV/c^{2}]");
  Arr.back().Sumw2();
  
  // mass13sq
  unsigned int sys13(helkin.addSubSystem({0}, {2}, {1}, {}));
  auto m13sq_limit = helkin.invMassBounds(sys13);
  double m13sq_min = m13sq_limit.first;
  double m13sq_max = m13sq_limit.second;
  Arr.push_back(
      TH1D("m13sq", "m_{13}^{2} [GeV/c^{2}]", NumBins, m13sq_min, m13sq_max));
  binWidth = (double)(m13sq_min - m13sq_max) / NumBins;
  sprintf(label, "Entries /%f.3", binWidth);
  Arr.back().GetYaxis()->SetTitle("# [" + TString(label) + "]");
  Arr.back().GetXaxis()->SetTitle("m_{13}^{2} [GeV/c^{2}]");
  Arr.back().Sumw2();
  
  // mass12sq
  unsigned int sys12(helkin.addSubSystem({0}, {1}, {2}, {}));
  auto m12sq_limit = helkin.invMassBounds(sys12);
  double m12sq_min = m12sq_limit.first;
  double m12sq_max = m12sq_limit.second;
  Arr.push_back(
      TH1D("m12sq", "m_{12}^{2} [GeV/c^{2}]", NumBins, m12sq_min, m12sq_max));
  binWidth = (double)(m12sq_min - m12sq_max) / NumBins;
  sprintf(label, "Entries /%f.3", binWidth);
  Arr.back().GetYaxis()->SetTitle("# [" + TString(label) + "]");
  Arr.back().GetXaxis()->SetTitle("m_{12}^{2} [GeV/c^{2}]");
  Arr.back().Sumw2();

  Arr2D.push_back(TH2D(TString(name + "_m23sqm13sq"), TString(title), NumBins,
                       m23sq_min, m23sq_max, NumBins, m13sq_min, m13sq_max));
  Arr2D.push_back(TH2D(TString(name + "_m23sqm12sq"), TString(title), NumBins,
                       m23sq_min, m23sq_max, NumBins, m12sq_min, m12sq_max));
  Arr2D.push_back(TH2D(TString(name + "_m12sqm13sq"), TString(title), NumBins,
                       m12sq_min, m12sq_max, NumBins, m13sq_min, m13sq_max));
  Arr2D.push_back(TH2D(TString(name + "_m23sqCosTheta"), TString(title),
                       NumBins, m23sq_min, m23sq_max, NumBins, -1, 1));

  Arr2D.at(0).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
  Arr2D.at(0).GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
  Arr2D.at(1).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}");
  Arr2D.at(1).GetYaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
  Arr2D.at(2).GetXaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
  Arr2D.at(2).GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
  Arr2D.at(3).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
  Arr2D.at(3).GetYaxis()->SetTitle("#cos(#Theta)_{KK}");

  auto itr = Arr2D.begin();
  for (; itr != Arr2D.end(); ++itr) {
    (*itr).GetXaxis()->SetNdivisions(508);
    (*itr).GetZaxis()->SetTitle("Entries");
  }

  setColor(color);
  return;
}
void DalitzHisto::fill(const HelicityKinematics &helkin,
                       const Data::DataSet &sample) {
  if(!sample.Data.size())
    throw std::runtime_error("DalitzHist::fill() | Empty data sample.");
  std::vector<double> w(sample.Weights.size(), 1.0);
  fill(helkin, sample, w);
}

void DalitzHisto::fill(const HelicityKinematics &helkin,
                       const Data::DataSet &sample, std::vector<double> w) {
  if(!sample.Data.size())
    throw std::runtime_error("DalitzHist::fill() | Empty data sample.");
  if(w.size() != sample.Weights.size())
    throw std::runtime_error("DalitzHist::fill() | Vector of weights and "
                             "sample do not have the same length");
  
  double weight =
      std::accumulate(sample.Weights.begin(), sample.Weights.end(), 0.0);

  Integral += weight;

  unsigned int sysId23 =
      helkin.getDataID(Physics::SubSystem({{1}, {2}}, {0}, {}));
  unsigned int sysId13 =
      helkin.getDataID(Physics::SubSystem({{0}, {2}}, {1}, {}));
  unsigned int sysId12 =
      helkin.getDataID(Physics::SubSystem({{0}, {1}}, {2}, {}));

  auto m23sq = sample.Data.at(3 * sysId23);
  auto cos23 = sample.Data.at(3 * sysId23 + 1);
  auto m13sq = sample.Data.at(3 * sysId13);
  auto m12sq = sample.Data.at(3 * sysId12);

  for (size_t i = 0; i < w.size(); ++i) {
    auto ww = sample.Weights.at(i)*w.at(i);
    Arr.at(0).Fill(m23sq.at(i), ww);
    Arr.at(1).Fill(m13sq.at(i), ww);
    Arr.at(2).Fill(m12sq.at(i), ww);

    Arr2D.at(0).Fill(m23sq.at(i), m13sq.at(i), ww);
    Arr2D.at(1).Fill(m23sq.at(i), m12sq.at(i), ww);
    Arr2D.at(2).Fill(m12sq.at(i), m13sq.at(i), ww);
    Arr2D.at(3).Fill(m23sq.at(i), cos23.at(i), ww);
    BranchPoint = {m23sq.at(i), m13sq.at(i), m12sq.at(i)};
    BranchWeight = ww;
    BranchEff = 1.0;
    Tree->Fill();
  }
}

void DalitzHisto::fill(const HelicityKinematics &helkin, const DataPoint &point,
                       double w) {

  double weight = point.Weight * w; // use event weights?

  Integral += weight;

  unsigned int sysId23 =
      helkin.getDataID(Physics::SubSystem({{1}, {2}}, {0}, {}));
  unsigned int sysId13 =
      helkin.getDataID(Physics::SubSystem({{0}, {2}}, {1}, {}));
  unsigned int sysId12 =
      helkin.getDataID(Physics::SubSystem({{0}, {1}}, {2}, {}));

  double m23sq = point.KinematicVariableList[3 * sysId23];
  double cos23 = point.KinematicVariableList[3 * sysId23 + 1];
  double m13sq = point.KinematicVariableList[3 * sysId13];
  //  double cos13 = point.getVal(3*sysId13+1);
  double m12sq = point.KinematicVariableList[3 * sysId12];
  //  double cos12 = point.getVal(3*sysId12+1);

  Arr.at(0).Fill(m23sq, weight);
  Arr.at(1).Fill(m13sq, weight);
  Arr.at(2).Fill(m12sq, weight);

  Arr2D.at(0).Fill(m23sq, m13sq, weight);
  Arr2D.at(1).Fill(m23sq, m12sq, weight);
  Arr2D.at(2).Fill(m12sq, m13sq, weight);
  Arr2D.at(3).Fill(m23sq, cos23, weight);
}

void DalitzHisto::setStats(bool b) {
  auto n = Arr.size();
  for (unsigned int i = 0; i < n; ++i) {
    Arr.at(i).SetStats(b);
  }
  auto n2 = Arr2D.size();
  for (unsigned int i = 0; i < n2; ++i) {
    Arr2D.at(i).SetStats(b);
  }
}

void DalitzHisto::scale(double w) {
  auto n = Arr.size();
  for (unsigned int i = 0; i < n; ++i) {
    Arr.at(i).Scale(w);
  }
  auto n2 = Arr2D.size();
  for (unsigned int i = 0; i < n2; ++i) {
    Arr2D.at(i).Scale(w);
  }
}

void DalitzHisto::setColor(Color_t color) {
  auto n = Arr.size();
  for (unsigned int i = 0; i < n; ++i) {
    Arr.at(i).SetLineColor(color);
    Arr.at(i).SetMarkerColor(color);
  }
}

TH1D *DalitzHisto::getHistogram(unsigned int num) { return &Arr.at(num); }

TH2D *DalitzHisto::getHistogram2D(unsigned int num) { return &Arr2D.at(num); }

void DalitzHisto::write() {
  Tree->Write(TString(Name) + "Tree");
  gDirectory->mkdir(TString(Name) + "_hist");
  gDirectory->cd(TString(Name) + "_hist");
  auto n = Arr.size();
  for (unsigned int i = 0; i < n; ++i) {
    Arr.at(i).Write();
  }
  auto n2 = Arr2D.size();
  for (unsigned int i = 0; i < n2; ++i) {
    Arr2D.at(i).Write();
  }
  gDirectory->cd("../");
}

} // namespace Plotting
} // namespace Tools
} // namespace ComPWA
