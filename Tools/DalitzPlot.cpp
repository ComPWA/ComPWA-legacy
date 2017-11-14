// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <stdio.h>
#include <numeric>

#include "Math/ProbFuncMathCore.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"

#include "Core/DataPoint.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Physics/Resonance.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

#include "Tools/DalitzPlot.hpp"

#include "Tools/HistTools.hpp"

using namespace ComPWA;
using namespace ComPWA::Physics::HelicityFormalism;
using namespace ComPWA::Tools;

void phspContour(unsigned int xsys, unsigned int ysys, unsigned int n,
                 double *xcoord, double *ycoord) {

  unsigned int num = n;
  if (num % 2 != 0) {
    num -= 1;
    BOOST_LOG_TRIVIAL(info)
        << "DalitzKinematics::phspContour() | "
           "Setting size to a even number. Assure that the size of "
           "your arrays is "
        << num * 2 + 1 << "!";
  }

  //  std::pair<double,double> xlimits = GetMinMax(xsys);
  //  double binw = (xlimits.second-xlimits.first)/(double)(num);
  //  double ymin,ymax,x;
  //  unsigned int i=0;
  //
  //  for (; i<num; i++)
  //    {
  //    x = i*binw + xlimits.first;
  //    while(x<=xlimits.first) { x+=binw/100; }
  //    while(x>=xlimits.second) { x-=binw/100; }
  //    auto lim = GetMinMaxLocal(ysys,xsys,x);
  //    //      ymin = invMassMin(ysys,xsys,x);
  //    //      ymax = invMassMax(ysys,xsys,x);
  //
  //    xcoord[i]=x; ycoord[i]=lim.first;
  //    xcoord[num*2-i]=x; ycoord[num*2-i]=lim.second;
  //    }
  //  //adding last datapoint
  //  x = i*binw + xlimits.first;
  //  while(x<=xlimits.first) { x+=binw/100; }
  //  while(x>=xlimits.second) { x-=binw/100; }
  //
  //  xcoord[i]=x; ycoord[i]=GetMinMaxLocal(ysys,xsys,x).first;
  return;
}

DalitzPlot::DalitzPlot(std::shared_ptr<Kinematics> kin, std::string name,
                       int bins)
    : _name(name), kin_(kin), _isFilled(0), _bins(bins), _globalScale(1.0),
      _correctForEfficiency(false),
      h_weights("h_weights", "h_weights", bins, 0, 1.01),
      dataDiagrams(kin, "data", "Data", bins),
      phspDiagrams(kin, "phsp", "Phase-space", bins),
      fitHitMissDiagrams(kin, "fitHitMiss", "HitMiss", bins) {
  gStyle->SetOptStat(10); // entries only
  //	gStyle->SetOptStat(1000001); //name and integral
  gStyle->SetOptTitle(0);

  // Full intensity blue
  // Phase-space green
  phspDiagrams.setColor(kGreen);

  phspDiagrams.SetStats(0);

  //=== generate contour
  double xpoints[4001], ypoints[4001];
  phspContour(0, 1, 2000, xpoints, ypoints);
  m23m13_contour = TGraph(4001, xpoints, ypoints);
  m23m13_contour.SetMarkerStyle(1);
  m23m13_contour.SetLineColor(kRed);
  m23m13_contour.SetMarkerColor(kRed);
  m23m13_contour.SetMarkerSize(0.0);
  m23m13_contour.SetTitle("phspContour");
  m23m13_contour.SetFillColor(kWhite);
  phspContour(0, 2, 2000, xpoints, ypoints);
  m23m12_contour = TGraph(4001, xpoints, ypoints);
  m23m12_contour.SetMarkerStyle(1);
  m23m12_contour.SetLineColor(kRed);
  m23m12_contour.SetMarkerColor(kRed);
  m23m12_contour.SetMarkerSize(0.0);
  m23m12_contour.SetTitle("phspContour");
  m23m12_contour.SetFillColor(kWhite);
  phspContour(2, 1, 2000, xpoints, ypoints);
  m12m13_contour = TGraph(4001, xpoints, ypoints);
  m12m13_contour.SetMarkerStyle(1);
  m12m13_contour.SetLineColor(kRed);
  m12m13_contour.SetMarkerColor(kRed);
  m12m13_contour.SetMarkerSize(0.0);
  m12m13_contour.SetTitle("phspContour");
  m12m13_contour.SetFillColor(kWhite);
}

DalitzPlot::~DalitzPlot() {}

void DalitzPlot::SetFitAmp(std::shared_ptr<ComPWA::AmpIntensity> intens,
                           std::string title, Color_t color) {
  _plotComponents.clear();
  _plotComponents.push_back(intens);
  _plotHistograms.push_back(
      DalitzHisto(kin_, intens->Name(), title, _bins, color));
  _plotHistograms.back().SetStats(0);
  _plotLegend.push_back("Fit");
  _isFilled = 0;
}

void DalitzPlot::Fill(std::shared_ptr<Kinematics> kin) {
  // TODO: reset diagrams here

  //===== Fill data histograms
  if (s_data) {
    for (unsigned int i = 0; i < s_data->GetNEvents(); i++) { // loop over data
      Event event(s_data->GetEvent(i));

      double eff = 1.0;
      if (_correctForEfficiency)
        eff = event.GetEfficiency();
      if (eff == 0.0) {
        LOG(error) << "DalitzPlot::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      double evWeight = event.GetWeight();

      dataDiagrams.Fill(kin, event, evWeight * 1 / eff);
      h_weights.Fill(evWeight * 1 / eff);
    }
    _globalScale = dataDiagrams.GetIntegral();
  }

  //===== Plot amplitude
  if (s_phsp) {
    LOG(info)
        << "PlotData::plot | Plotting phase space sample and intensity...";

    double weightsSum = 0.0;

    // Loop over all events in phase space sample
    progressBar bar(s_phsp->GetNEvents());
    for (unsigned int i = 0; i < s_phsp->GetNEvents();
         i++) { // loop over phsp MC
      bar.nextEvent();
      Event event(s_phsp->GetEvent(i));
      double evWeight = event.GetWeight();
      double eff = 1.0;
      if (_correctForEfficiency)
        eff = event.GetEfficiency();
      if (eff == 0.0) {
        LOG(error) << "DalitzPlot::Fill() | Loop over "
                      "phsp sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }

      double evBase = evWeight / eff;
      weightsSum += evBase;

      /* Construct DataPoint from Event to check if dataPoint is within the
       * phase space boundaries */
      dataPoint point;
      try {
        kin->EventToDataPoint(event, point);
      } catch (BeyondPhsp &ex) { // event outside phase, remove
        continue;
      }

      // Fill diagrams with pure phase space events
      phspDiagrams.Fill(kin, event, evBase); // scale phsp to data size

      // Loop over all components that we want to plot
      for (int t = 0; t < _plotHistograms.size(); t++)
        _plotHistograms.at(t).Fill(
            kin, event, _plotComponents.at(t)->Intensity(point) * evBase);
    }

    // Scale histograms to match data sample
    phspDiagrams.Scale(_globalScale / phspDiagrams.GetIntegral());
    double scale = _globalScale / _plotHistograms.at(0).GetIntegral();
    _plotHistograms.at(0).Scale(scale);

    for (int t = 1; t < _plotHistograms.size(); t++) {
      _plotHistograms.at(t).Scale(scale);
    }
  }

  //===== Plot hit&miss data
  if (s_hitMiss) {
    for (unsigned int i = 0; i < s_hitMiss->GetNEvents();
         i++) { // loop over data
      Event event(s_hitMiss->GetEvent(i));
      double eff = 1.0;
      if (_correctForEfficiency)
        eff = event.GetEfficiency();
      if (eff == 0.0) {
        LOG(error) << "DalitzPlot::Fill() | Loop over "
                      "Hit&Miss sample: An event with zero efficiency was "
                      "found! This should not happen! We skip it!";
        continue;
      }
      double evWeight = event.GetWeight();

      fitHitMissDiagrams.Fill(kin, event, evWeight * 1 / eff);
    }
  }

  _isFilled = 1;
}

void DalitzPlot::Plot() {
  if (!_isFilled)
    Fill(kin_);

  //----- plotting 2D dalitz distributions -----
  TCanvas *c1 = new TCanvas("dalitz", "dalitz", 50, 50, 1600, 1600);
  c1->Divide(2, 2);
  c1->cd(1);
  dataDiagrams.getHistogram2D(0)->Draw("COLZ");
  m23m13_contour.Draw("P");
  dataDiagrams.getHistogram2D(0)->SetLineWidth(1);
  c1->cd(2);
  fitHitMissDiagrams.getHistogram2D(0)->Draw("COLZ");
  m23m13_contour.Draw("P");
  c1->cd(3);
  _plotHistograms.at(0).getHistogram2D(0)->Draw("COLZ");
  m23m13_contour.Draw("P");

  //----- plotting invariant mass distributions -----
  TCanvas *c2 = new TCanvas("invmass", "invmass", 50, 50, 2400, 800);
  c2->Divide(3, 1);
  c2->cd(1);
  CreateHist(0); // Plotting mKKsq
  c2->cd(2);
  CreateHist(1); // Plotting mKSK+sq
  c2->cd(3);
  CreateHist(2); // Plotting mKSK+sq
  c2->cd(3);
  TLegend *leg = new TLegend(0.15, 0.6, 0.50, 0.85);
  leg->AddEntry(dataDiagrams.getHistogram(2), "Data");
  leg->AddEntry(_plotHistograms.at(0).getHistogram(2), "Model");
  for (int i = 1; i < _plotComponents.size(); i++)
    leg->AddEntry(_plotHistograms.at(i).getHistogram(2),
                  TString(_plotLegend.at(i)));
  leg->SetFillStyle(0);
  leg->Draw(); // Plot legend

  //----- plotting signal amplitude contributions -----
  TCanvas *c5 =
      new TCanvas("signalInvmass", "signalInvmass", 50, 50, 2400, 800);
  c5->Divide(3, 1);
  c5->cd(1);
  CreateHist2(0); // Plotting mKKsq
  c5->cd(2);
  CreateHist2(1); // Plotting mKSK+sq
  c5->cd(3);
  CreateHist2(2); // Plotting mKSK+sq

  //----- Helicity angle distributions -----
  //  TCanvas *c3 =
  //      new TCanvas("helicityAngle", "helicity angle", 50, 50, 2400, 1600);
  //  c3->Divide(3, 2);
  //  c3->cd(1);
  //  CreateHist(3);
  //  c3->cd(2);
  //  CreateHist(4);
  //  c3->cd(3);
  //  CreateHist(5);
  //  c3->cd(4);
  //  CreateHist(6);
  //  c3->cd(5);
  //  CreateHist(7);
  //  c3->cd(6);
  //  CreateHist(8);

  //----- Weights distributions -----
  TCanvas *c4 = new TCanvas("weights", "weights", 50, 50, 2400, 800);
  c4->cd();
  h_weights.Draw();

  //----- Write to TFile -----
  TFile *tf2 = new TFile(_name + ".root", "recreate");
  if (tf2->IsZombie()) {
    std::cout << "Error opening output file" << std::endl;
    exit(-1);
  }
  m12m13_contour.Write("m12m13_contour", TObject::kOverwrite, 0);
  m23m12_contour.Write("m23m12_contour", TObject::kOverwrite, 0);
  m23m13_contour.Write("m23m13_contour", TObject::kOverwrite, 0);
  c1->Write("dalitz", TObject::kOverwrite, 0);
  c2->Write("invmass", TObject::kOverwrite, 0);
  c5->Write("signalInvmass", TObject::kOverwrite, 0);
  //  c3->Write("helicityAngle", TObject::kOverwrite, 0);

  // Save data trees and histograms
  tf2->mkdir("hist");
  tf2->cd("hist");
  dataDiagrams.Write();
  phspDiagrams.Write();
  fitHitMissDiagrams.Write();

  for (int t = 0; t < _plotHistograms.size(); t++)
    _plotHistograms.at(t).Write();

  // Write some canvas to single files
  c2->Print(_name + "-invmass.root");
  c2->Print(_name + "-invmass.pdf");

  tf2->Close();

  return;
}

void DalitzPlot::CreateHist(unsigned int id) {
  TPad *pad;
  std::vector<TH1D *> v;
  std::vector<TString> options;
  if (s_data) {
    v.push_back(dataDiagrams.getHistogram(id));
    options.push_back("E1");
  }
  for (int t = 0; t < _plotHistograms.size(); t++) {
    v.push_back(_plotHistograms.at(t).getHistogram(id));
    options.push_back("Sames,Hist");
  }

  pad = drawPull(v, options);
}

void DalitzPlot::CreateHist2(unsigned int id) {
  TPad *pad;
  std::vector<TH1D *> v;
  std::vector<TString> options;

  v.push_back(_plotHistograms.at(0).getHistogram(id));
  options.push_back("Hist");
  for (int t = 1; t < _plotHistograms.size(); t++) {
    v.push_back(_plotHistograms.at(t).getHistogram(id));
    options.push_back("Sames,Hist");
  }

  pad = drawHist(v, options);
}

//===================== DalitzHisto =====================
DalitzHisto::DalitzHisto(std::shared_ptr<Kinematics> kin, std::string name,
                         std::string title, unsigned int bins, Color_t color)
    : _name(name), _title(title), _nBins(bins), _integral(0.0), _color(color) {

  // we have to explicitly cast to HelicityKinematics in order to get
  // the invariant mass boundaries
  auto helkin = std::dynamic_pointer_cast<HelicityKinematics>(kin);

  // Initialize TTree
  _tree = std::unique_ptr<TTree>(new TTree(TString(_name), TString(_title)));

  // Adding branches to TTree
  _tree->Branch(TString(_name), &t_point);
  _tree->Branch("efficiency", &t_eff, "eff/D");
  _tree->Branch("weight", &t_weight, "weight/D");

  char label[60];

  // mass23sq
  SubSystem sys23({0}, {1}, {2});
  auto m23sq_limit = helkin->GetInvMassBounds(sys23);
  double m23sq_min = m23sq_limit.first;
  double m23sq_max = m23sq_limit.second;

  _arr.push_back(
      TH1D("m23sq", "m_{23}^{2} [GeV/c^{2}]", _nBins, m23sq_min, m23sq_max));
  double binWidth = (double)(m23sq_min - m23sq_max) / _nBins;
  sprintf(label, "Entries /%f.3", binWidth);
  _arr.back().GetYaxis()->SetTitle("# [" + TString(label) + "]");
  _arr.back().GetXaxis()->SetTitle("m_{23}^{2} [GeV/c^{2}]");
  _arr.back().Sumw2();
  // mass13sq
  SubSystem sys13({1}, {0}, {2});
  auto m13sq_limit = helkin->GetInvMassBounds(sys13);
  double m13sq_min = m13sq_limit.first;
  double m13sq_max = m13sq_limit.second;

  _arr.push_back(
      TH1D("m13sq", "m_{13}^{2} [GeV/c^{2}]", _nBins, m13sq_min, m13sq_max));
  binWidth = (double)(m13sq_min - m13sq_max) / _nBins;
  sprintf(label, "Entries /%f.3", binWidth);
  _arr.back().GetYaxis()->SetTitle("# [" + TString(label) + "]");
  _arr.back().GetXaxis()->SetTitle("m_{13}^{2} [GeV/c^{2}]");
  _arr.back().Sumw2();
  // mass12sq
  SubSystem sys12({0}, {0}, {1});
  auto m12sq_limit = helkin->GetInvMassBounds(sys12);
  double m12sq_min = m12sq_limit.first;
  double m12sq_max = m12sq_limit.second;

  _arr.push_back(
      TH1D("m12sq", "m_{12}^{2} [GeV/c^{2}]", _nBins, m12sq_min, m12sq_max));
  binWidth = (double)(m12sq_min - m12sq_max) / _nBins;
  sprintf(label, "Entries /%f.3", binWidth);
  _arr.back().GetYaxis()->SetTitle("# [" + TString(label) + "]");
  _arr.back().GetXaxis()->SetTitle("m_{12}^{2} [GeV/c^{2}]");
  _arr.back().Sumw2();

  _arr2D.push_back(TH2D(TString(name + "_m23sqm13sq"), TString(title), _nBins,
                        m23sq_min, m23sq_max, _nBins, m13sq_min, m13sq_max));
  _arr2D.push_back(TH2D(TString(name + "_m23sqm12sq"), TString(title), _nBins,
                        m23sq_min, m23sq_max, _nBins, m12sq_min, m12sq_max));
  _arr2D.push_back(TH2D(TString(name + "_m12sqm13sq"), TString(title), _nBins,
                        m12sq_min, m12sq_max, _nBins, m13sq_min, m13sq_max));
  _arr2D.push_back(TH2D(TString(name + "_m23sqCosTheta"), TString(title),
                        _nBins, m23sq_min, m23sq_max, _nBins, -1, 1));

  _arr2D.at(0).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
  _arr2D.at(0).GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
  _arr2D.at(1).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}");
  _arr2D.at(1).GetYaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
  _arr2D.at(2).GetXaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
  _arr2D.at(2).GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
  _arr2D.at(3).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
  _arr2D.at(3).GetYaxis()->SetTitle("#cos(#Theta)_{KK}");

  auto itr = _arr2D.begin();
  for (; itr != _arr2D.end(); ++itr) {
    (*itr).GetXaxis()->SetNdivisions(508);
    (*itr).GetZaxis()->SetTitle("Entries");
  }

  setColor(_color);
  return;
}

void DalitzHisto::Fill(std::shared_ptr<Kinematics> kin, Event &event,
                       double w) {

  double weight = event.GetWeight() * w; // use event weights?

  _integral += weight;

  int sysId23 = kin->GetDataID(SubSystem({0}, {1}, {2}));
  int sysId13 = kin->GetDataID(SubSystem({1}, {0}, {2}));
  int sysId12 = kin->GetDataID(SubSystem({2}, {0}, {1}));

  dataPoint point;
  try {
    kin->EventToDataPoint(event, point);
  } catch (std::exception &ex) {
    return;
  }

  double m23sq = point.GetValue(3 * sysId23);
  double cos23 = point.GetValue(3 * sysId23 + 1);
  double m13sq = point.GetValue(3 * sysId13);
  //	double cos13 = point.getVal(3*sysId13+1);
  double m12sq = point.GetValue(3 * sysId12);
  //	double cos12 = point.getVal(3*sysId12+1);

  _arr.at(0).Fill(m23sq, weight);
  _arr.at(1).Fill(m13sq, weight);
  _arr.at(2).Fill(m12sq, weight);

  _arr2D.at(0).Fill(m23sq, m13sq, weight);
  _arr2D.at(1).Fill(m23sq, m12sq, weight);
  _arr2D.at(2).Fill(m12sq, m13sq, weight);
  _arr2D.at(3).Fill(m23sq, cos23, weight);
}

void DalitzHisto::SetStats(bool b) {
  auto n = _arr.size();
  for (int i = 0; i < n; ++i) {
    _arr.at(i).SetStats(b);
  }
  auto n2 = _arr2D.size();
  for (int i = 0; i < n2; ++i) {
    _arr2D.at(i).SetStats(b);
  }
}

void DalitzHisto::Scale(double w) {
  auto n = _arr.size();
  for (int i = 0; i < n; ++i) {
    _arr.at(i).Scale(w);
  }
  auto n2 = _arr2D.size();
  for (int i = 0; i < n2; ++i) {
    _arr2D.at(i).Scale(w);
  }
}

void DalitzHisto::setColor(Color_t color) {
  auto n = _arr.size();
  for (int i = 0; i < n; ++i) {
    _arr.at(i).SetLineColor(color);
    _arr.at(i).SetMarkerColor(color);
  }
}

TH1D *DalitzHisto::getHistogram(unsigned int num) { return &_arr.at(num); }

TH2D *DalitzHisto::getHistogram2D(unsigned int num) { return &_arr2D.at(num); }

void DalitzHisto::Write() {
  _tree->Write(TString(_name) + "_tree");
  gDirectory->mkdir(TString(_name) + "_hist");
  gDirectory->cd(TString(_name) + "_hist");
  auto n = _arr.size();
  for (int i = 0; i < n; ++i) {
    _arr.at(i).Write();
  }
  auto n2 = _arr2D.size();
  for (int i = 0; i < n2; ++i) {
    _arr2D.at(i).Write();
  }
  gDirectory->cd("../");
}
