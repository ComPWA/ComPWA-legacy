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

#include "PlotData.hpp"

#include "Tools.hpp"

using namespace ComPWA;
using namespace ComPWA::Physics::HelicityFormalism;

void phspContour(unsigned int xsys,unsigned int ysys,
                                   unsigned int n, double* xcoord, double* ycoord)
{
  
  unsigned int num=n;
  if(num%2!=0) {
    num-=1;
    BOOST_LOG_TRIVIAL(info)<<"DalitzKinematics::phspContour() | "
    "Setting size to a even number. Assure that the size of "
    "your arrays is "<<num*2+1<<"!";
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

plotData::plotData(std::string name, int bins)
    : _name(name), _isFilled(0), _bins(bins), _globalScale(1.0),
      h_weights("h_weights", "h_weights", bins, 0, 1.01),
      dataDiagrams("data", "Data", bins),
      phspDiagrams("phsp", "Phase-space", bins),
      fitDiagrams("fit", "Model", bins),
      fitHitMissDiagrams("fitHitMiss", "HitMiss", bins) {
  gStyle->SetOptStat(10); // entries only
  //	gStyle->SetOptStat(1000001); //name and integral
  gStyle->SetOptTitle(0);

  // Full intensity blue
  fitDiagrams.setColor(kBlue - 4);
  // Phase-space green
  phspDiagrams.setColor(kGreen);

  fitDiagrams.SetStats(0);
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

plotData::~plotData() {}

void plotData::SetFitAmp(std::shared_ptr<AmpIntensity> intens) {
  _intens = intens;
  _isFilled = 0;
}

void plotData::Fill() {
  // TODO: reset diagrams here

  //===== Fill data histograms
  if (s_data) {
    for (unsigned int i = 0; i < s_data->getNEvents(); i++) { // loop over data
      Event event(s_data->getEvent(i));

      double eff = 1.0;
      if (_correctForEfficiency)
        eff = event.getEfficiency();
      if (eff == 0.0) {
        LOG(error) << "plotData::Fill() | Loop over "
                      "data sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      double evWeight = event.getWeight();

      dataDiagrams.Fill(event, evWeight * 1 / eff);
      h_weights.Fill(evWeight * 1 / eff);
    }
    _globalScale = dataDiagrams.GetIntegral();
  }

  //===== Plot amplitude
  LOG(info) << "PlotData::plot | Evaluating amplitude...";

  if (_intens && s_phsp) {

    /* Loop over all events in phase space sample */
    progressBar bar(s_phsp->getNEvents());
    for (unsigned int i = 0; i < s_phsp->getNEvents();
         i++) { // loop over phsp MC
      bar.nextEvent();
      Event event(s_phsp->getEvent(i));
      double eff = 1.0;
      if (_correctForEfficiency)
        eff = event.getEfficiency();
      if (eff == 0.0) {
        LOG(error) << "plotData::Fill() | Loop over "
                      "phsp sample: An event with zero efficiency was found! "
                      "This should not happen! We skip it!";
        continue;
      }
      //double evWeight = event.getWeight();

      /* Construct DataPoint from Event to check if dataPoint is within the
       * phase space boundaries */
      dataPoint point;
      try {
        point = dataPoint(event);
      } catch (BeyondPhsp &ex) { // event outside phase, remove
        continue;
      }

      // Fill diagrams with pure phase space events
      phspDiagrams.Fill(event, 1 / eff); // scale phsp to data size

//      for (int t = 0; t < signalComponents.size(); ++t) {
//        resonanceItr res = _ampVec.at(0)->GetResonanceItrFirst();
//        for (int j = 0; j < t; ++j)
//          res++;
//
//        // skip CP partner of resonance
//        if ((*res)->GetName().find("_CP") != std::string::npos)
//          continue;
//
//        std::complex<double> val = (*res)->Evaluate(point);
//
//        try { // trying to find a CP partner and add it
//          std::shared_ptr<Resonance> cpRes;
//          std::shared_ptr<ComPWA::Physics::AmplitudeSum::AmpSumIntensity>
//              tmpAmp = std::dynamic_pointer_cast<
//                  ComPWA::Physics::AmplitudeSum::AmpSumIntensity>(
//                  _ampVec.at(0));
//          cpRes = tmpAmp->GetResonance((*res)->GetName() + "_CP");
//          val += cpRes->Evaluate(point);
//        } catch (std::exception &ex) {
//        }
//        signalComponents.at(t).Fill(point, std::norm(val) / eff);
//      }

      /* Loop over all resonances of first amplitude. This is supposed to be our
       * signal intensity */
//      std::complex<double> tmp_intens2(0, 0);
//      auto it = _ampVec.at(0)->GetResonanceItrFirst();
//      for (; it != _ampVec.at(0)->GetResonanceItrLast(); ++it) {
//        if ((*it)->GetName().find("_CP") != std::string::npos)
//          continue;
//        tmp_intens2 += (*it)->Evaluate(point);
//      }
//      ampHistos.at(0).Fill(point, std::norm(tmp_intens2) / eff);

      /* Loop over all amplitudes. This is supposed to be our total intensity
       * (usually signal+background) */
      double intens = 0;
      intens += _intens->Intensity(point);
      fitDiagrams.Fill(event, intens / eff);
    }

    // Scale histograms to match data sample
    fitDiagrams.Scale(_globalScale / fitDiagrams.GetIntegral());
    phspDiagrams.Scale(_globalScale / phspDiagrams.GetIntegral());

//    for (int t = 0; t < _ampVec.size(); ++t) {
//      double scale =
//          _globalScale / ampHistos.at(t).GetIntegral() * _fraction.at(t);
//      ampHistos.at(t).Scale(scale);
//    }
//    for (int t = 0; t < signalComponents.size(); ++t) {
//      double scale =
//          _globalScale / ampHistos.at(0).GetIntegral() * _fraction.at(0);
//      signalComponents.at(t).Scale(scale);
//    }
  }

  //===== Plot hit&miss data
  if (s_hitMiss) {
    for (unsigned int i = 0; i < s_hitMiss->getNEvents();
         i++) { // loop over data
      Event event(s_hitMiss->getEvent(i));
      double eff = 1.0;
      if (_correctForEfficiency)
        eff = event.getEfficiency();
      if (eff == 0.0) {
        LOG(error) << "plotData::Fill() | Loop over "
                      "Hit&Miss sample: An event with zero efficiency was "
                      "found! This should not happen! We skip it!";
        continue;
      }
      double evWeight = event.getWeight();

      fitHitMissDiagrams.Fill(event, evWeight * 1 / eff);
    }
  }

  _isFilled = 1;
}

void plotData::Plot() {
  if (!_isFilled)
    Fill();

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
  fitDiagrams.getHistogram2D(0)->Draw("COLZ");
  m23m13_contour.Draw("P");

  //----- plotting invariant mass distributions -----
  TCanvas *c2 = new TCanvas("invmass", "invmass", 50, 50, 2400, 800);
  c2->Divide(3, 1);
    c2->cd(1);
    CreateHist(0);// Plotting mKKsq
    c2->cd(2);
    CreateHist(1);// Plotting mKSK+sq
    c2->cd(3);
    CreateHist(2);// Plotting mKSK+sq
    c2->cd(3);
  TLegend *leg = new TLegend(0.15, 0.6, 0.50, 0.85);
  leg->AddEntry(dataDiagrams.getHistogram(2), "Data");
  leg->AddEntry(fitDiagrams.getHistogram(2), "Model");
  if (ampHistos.size())
    leg->AddEntry(ampHistos.back().getHistogram(2), "Background");
  leg->SetFillStyle(0);
  leg->Draw(); // Plot legend

  //----- plotting signal amplitude contributions -----
  TCanvas *c5 =
      new TCanvas("signalInvmass", "signalInvmass", 50, 50, 2400, 800);
  c5->Divide(3, 1);
    c5->cd(1);
    CreateHist2(0);// Plotting mKKsq
    c5->cd(2);
    CreateHist2(1);// Plotting mKSK+sq
  c5->cd(3);
  CreateHist2(2);// Plotting mKSK+sq

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
  fitDiagrams.Write();
  fitHitMissDiagrams.Write();
  auto itAmp = ampHistos.begin();
  for (; itAmp != ampHistos.end(); ++itAmp)
    (*itAmp).Write();

  itAmp = signalComponents.begin();
  for (; itAmp != signalComponents.end(); ++itAmp)
    (*itAmp).Write();

  // Write some canvas to single files
  c2->Print(_name + "-invmass.root");
  c2->Print(_name + "-invmass.pdf");

  tf2->Close();

  return;
}

void plotData::CreateHist(unsigned int id) {
  TPad *pad;
  std::vector<TH1D *> v;
  std::vector<TString> options;
  if (s_data) {
    v.push_back(dataDiagrams.getHistogram(id));
    options.push_back("E1");
  }
//  if (_ampVec.size()) {
//    v.push_back(fitDiagrams.getHistogram(id));
//    options.push_back("Sames,Hist");
//    for (unsigned int i = 0; i < plotComponent.size(); ++i) {
//      ampHistos.at(plotComponent.at(i).first)
//          .setColor(plotComponent.at(i).second);
//      v.push_back(ampHistos.at(plotComponent.at(i).first).getHistogram(id));
//      options.push_back("Sames,Hist");
//    }
//  }
  pad = drawPull(v, options);
}

void plotData::CreateHist2(unsigned int id) {
  TPad *pad;
  std::vector<TH1D *> v;
  std::vector<TString> options;
  //	v.push_back(ampHistos.at(0).getHistogram(id));
  v.push_back(fitDiagrams.getHistogram(id));
  options.push_back("Hist");

  if (signalComponents.size()) {
    for (unsigned int i = 0; i < signalComponents.size(); ++i) {
      v.push_back(signalComponents.at(i).getHistogram(id));
      options.push_back("Sames,Hist");
    }
  }
  pad = drawHist(v, options);
}

//===================== dalitzHisto =====================
dalitzHisto::dalitzHisto(std::string n, std::string t, unsigned int bins)
    : name(n), title(t), nBins(bins), _integral(0.0) {
  auto *kin = dynamic_cast<HelicityKinematics *>( Kinematics::Instance() );

  // Initialize TTree
  tree = std::unique_ptr<TTree>(new TTree(TString(name), TString(title)));

  // Adding branches to TTree
  tree->Branch(TString(name), &t_point);
  tree->Branch("efficiency", &t_eff, "eff/D");
  tree->Branch("weight", &t_weight, "weight/D");

  char label[60];

      //mass23sq
  SubSystem sys23({0}, {1}, {2});
  auto m23sq_limit = kin->GetInvMassBounds(sys23);
  double m23sq_min = m23sq_limit.first;
  double m23sq_max = m23sq_limit.second;

  arr.push_back(TH1D("m23sq", "m_{23}^{2} [GeV/c^{2}]", nBins,
                     m23sq_min, m23sq_max));
  double binWidth = (double)(m23sq_min - m23sq_max) / nBins;
  sprintf(label, "Entries /%f.3", binWidth);
  arr.back().GetYaxis()->SetTitle("# [" + TString(label) + "]");
  arr.back().GetXaxis()->SetTitle("m_{23}^{2} [GeV/c^{2}]");
  arr.back().Sumw2();
      //mass13sq
  SubSystem sys13({1}, {0}, {2});
  auto m13sq_limit = kin->GetInvMassBounds(sys13);
  double m13sq_min = m13sq_limit.first;
  double m13sq_max = m13sq_limit.second;

  arr.push_back(TH1D("m13sq", "m_{13}^{2} [GeV/c^{2}]", nBins,
                     m13sq_min, m13sq_max));
  binWidth = (double)(m13sq_min - m13sq_max) / nBins;
  sprintf(label, "Entries /%f.3", binWidth);
  arr.back().GetYaxis()->SetTitle("# [" + TString(label) + "]");
  arr.back().GetXaxis()->SetTitle("m_{13}^{2} [GeV/c^{2}]");
  arr.back().Sumw2();
        //mass12sq
  SubSystem sys12({0}, {0}, {1});
  auto m12sq_limit = kin->GetInvMassBounds(sys12);
  double m12sq_min = m12sq_limit.first;
  double m12sq_max = m12sq_limit.second;

  arr.push_back(TH1D("m12sq", "m_{12}^{2} [GeV/c^{2}]", nBins,
                     m12sq_min, m12sq_max));
  binWidth = (double)(m12sq_min - m12sq_max) / nBins;
  sprintf(label, "Entries /%f.3", binWidth);
  arr.back().GetYaxis()->SetTitle("# [" + TString(label) + "]");
  arr.back().GetXaxis()->SetTitle("m_{12}^{2} [GeV/c^{2}]");
  arr.back().Sumw2();
      
  arr2D.push_back(TH2D(TString(name + "_m23sqm13sq"), TString(title), nBins,
                       m23sq_min, m23sq_max, nBins, m13sq_min, m13sq_max));
  arr2D.push_back(TH2D(TString(name + "_m23sqm12sq"), TString(title), nBins,
                       m23sq_min, m23sq_max, nBins, m12sq_min, m12sq_max));
  arr2D.push_back(TH2D(TString(name + "_m12sqm13sq"), TString(title), nBins,
                       m12sq_min, m12sq_max, nBins, m13sq_min, m13sq_max));
  arr2D.push_back(TH2D(TString(name + "_m23sqCosTheta"), TString(title), nBins,
                       m23sq_min, m23sq_max, nBins, -1, 1));

  arr2D.at(0).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
  arr2D.at(0).GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
  arr2D.at(1).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}");
  arr2D.at(1).GetYaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
  arr2D.at(2).GetXaxis()->SetTitle("m_{K_{S}K^{-}}^{2} [GeV^{2}/c^{4}]");
  arr2D.at(2).GetYaxis()->SetTitle("m_{K_{S}K^{+}}^{2} [GeV^{2}/c^{4}]");
  arr2D.at(3).GetXaxis()->SetTitle("m_{KK}^{2} [GeV^{2}/c^{4}]");
  arr2D.at(3).GetYaxis()->SetTitle("#cos(#Theta)_{KK}");

  auto itr = arr2D.begin();
  for (; itr != arr2D.end(); ++itr) {
    (*itr).GetXaxis()->SetNdivisions(508);
    (*itr).GetZaxis()->SetTitle("Entries");
  }
  return;
}
void dalitzHisto::Fill(Event &event, double w) {

  double weight = event.getWeight() * w; // use event weights?

  _integral += weight;
  
  auto *kin = dynamic_cast<HelicityKinematics*>( Kinematics::Instance() );
  
  SubSystem sys23({0}, {1}, {2});
  SubSystem sys13({1}, {0}, {2});
  SubSystem sys12({2}, {0}, {1});
  
  dataPoint point;
  kin->EventToDataPoint(event, point, sys23);
  kin->EventToDataPoint(event, point, sys13);
  kin->EventToDataPoint(event, point, sys12);
  
  double m23sq = point.GetValue(0);
  double cos23 = point.GetValue(1);
  double m13sq = point.GetValue(3);
  //	double cos13 = point.getVal(4);
  double m12sq = point.GetValue(6);
  //	double cos12 = point.getVal(7);
  
  arr.at(0).Fill(m23sq, weight);
  arr.at(1).Fill(m13sq, weight);
  arr.at(2).Fill(m12sq, weight);
  
  arr2D.at(0).Fill(m23sq, m13sq, weight);
  arr2D.at(1).Fill(m23sq, m12sq, weight);
  arr2D.at(2).Fill(m12sq, m13sq, weight);
  arr2D.at(3).Fill(m23sq, cos23, weight);
}

void dalitzHisto::SetStats(bool b) {
  auto n = arr.size();
  for (int i = 0; i < n; ++i) {
    arr.at(i).SetStats(b);
  }
  auto n2 = arr2D.size();
  for (int i = 0; i < n2; ++i) {
    arr2D.at(i).SetStats(b);
  }
}

void dalitzHisto::Scale(double w) {
  auto n = arr.size();
  for (int i = 0; i < n; ++i) {
    arr.at(i).Scale(w);
  }
  auto n2 = arr2D.size();
  for (int i = 0; i < n2; ++i) {
    arr2D.at(i).Scale(w);
  }
}

void dalitzHisto::setColor(Color_t color) {
  auto n = arr.size();
  for (int i = 0; i < n; ++i) {
    arr.at(i).SetLineColor(color);
    arr.at(i).SetMarkerColor(color);
  }
}

TH1D *dalitzHisto::getHistogram(unsigned int num) { return &arr.at(num); }

TH2D *dalitzHisto::getHistogram2D(unsigned int num) { return &arr2D.at(num); }

TH2Poly *dalitzHisto::getTH2PolyPull(TH2Poly *hist1, TH2Poly *hist2,
                                     TString name) {
  if (hist1->GetBins()->GetEntries() != hist2->GetBins()->GetEntries()) {
    std::cout << "binning doesnt match" << std::endl;
    return 0;
  }
  TH2Poly *resHist = (TH2Poly *)hist1->Clone(name + hist1->GetName());
  resHist->Reset("");

  hist2->Scale(hist1->Integral() / hist2->Integral());
  double limit = 0;
  for (int bin = 1; bin <= hist1->GetBins()->GetEntries(); bin++) {
    double c1 = hist1->GetBinContent(bin);
    double c2 = hist2->GetBinContent(bin);
    double c1Err = hist1->GetBinError(bin);
    double c2Err = hist2->GetBinError(bin);
    if (c1 > 0 && c2 > 0) {
      double pull = (c1 - c2) / sqrt(c1Err * c1Err + c2Err * c2Err);
      resHist->SetBinContent(bin, pull);
      if (std::abs(pull) > limit)
        limit = std::abs(pull);
    } else
      resHist->SetBinContent(bin, -999); // empty bins are set to error value
    resHist->SetBinError(bin, 0);
  }
  limit = 4;
  limit += 1;
  resHist->GetZaxis()->SetRangeUser(
      -limit, limit); // symmetric range including all values
  resHist->GetZaxis()->SetTitle("deviation [#sigma]");
  resHist->GetZaxis()->SetTitleOffset(0.5);

  return resHist;
}

void dalitzHisto::Write() {
  tree->Write(TString(name) + "_tree");
  gDirectory->mkdir(TString(name) + "_hist");
  gDirectory->cd(TString(name) + "_hist");
  auto n = arr.size();
  for (int i = 0; i < n; ++i) {
    arr.at(i).Write();
  }
  auto n2 = arr2D.size();
  for (int i = 0; i < n2; ++i) {
    arr2D.at(i).Write();
  }
  gDirectory->cd("../");
}
