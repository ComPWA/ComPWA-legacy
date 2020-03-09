// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Collection of useful routines to plot histograms,
/// e.g. creating (normalized) residual plots
///

#ifndef TOOLS_H
#define TOOLS_H

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TKDTreeBinning.h>
#include <TLatex.h>
#include <TMath.h>
#include <TObject.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TRandom3.h>
#include <TStyle.h>

#include "Core/Logging.hpp"
inline TH1 *getPull(TH1 *Histogram1, TH1 *Histogram2, TString Name = "pull_") {
  if (Histogram1->GetNbinsX() != Histogram2->GetNbinsX() ||
      Histogram1->GetNbinsY() != Histogram2->GetNbinsY() ||
      Histogram1->GetNbinsZ() != Histogram2->GetNbinsZ()) {
    std::cout << "binning doesn't match" << std::endl;
    return 0;
  }
  TH1 *PullHistogram = (TH1 *)Histogram1->Clone(Name + Histogram1->GetName());
  PullHistogram->Reset();
  double PullLimit = 0;
  for (int i = 1; i <= Histogram1->GetNbinsX(); ++i) {
    for (int j = 1; j <= Histogram1->GetNbinsY(); ++j) {
      for (int k = 1; k <= Histogram1->GetNbinsZ(); ++k) {
        unsigned int bin = Histogram1->GetBin(i, j, k);
        double c1 = Histogram1->GetBinContent(bin);
        double c2 = Histogram2->GetBinContent(bin);
        double c1Err = Histogram1->GetBinError(bin);
        double c2Err = Histogram2->GetBinError(bin);
        double Pull;
        if (c1 > 0 && c2 > 0) {
          Pull = (c1 - c2) / sqrt(c1Err * c1Err + c2Err * c2Err);
          PullHistogram->SetBinContent(bin, Pull);
          if (std::fabs(Pull) > PullLimit)
            PullLimit = std::abs(Pull);
        } else
          PullHistogram->SetBinContent(
              bin,
              -999); // empty bins are set to error value
        PullHistogram->SetBinError(bin, .0001);
      }
    }
  }
  PullLimit = 4; // set fix limits
  PullLimit += 1;
  // X
  PullHistogram->GetXaxis()->SetTitleSize(.14);
  PullHistogram->GetXaxis()->SetTitleOffset(.93);
  PullHistogram->GetXaxis()->SetLabelSize(.12);
  if (PullHistogram->GetDimension() == 1) {
    PullHistogram->GetYaxis()->SetRangeUser(-PullLimit, PullLimit);
    PullHistogram->GetYaxis()->SetTitle("deviation [#sigma]");
    PullHistogram->GetYaxis()->SetTitleOffset(0.36);
    PullHistogram->GetYaxis()->SetTitleSize(0.12);
    PullHistogram->GetYaxis()->SetLabelSize(0.12);
    PullHistogram->GetYaxis()->SetNdivisions(504);
    PullHistogram->GetYaxis()->CenterTitle();
  }
  if (PullHistogram->GetDimension() ==
      2) { // symmetric range including all values
    PullHistogram->GetZaxis()->SetRangeUser(-PullLimit, PullLimit);
    PullHistogram->GetZaxis()->SetTitle("deviation [#sigma]");
    PullHistogram->GetYaxis()->SetTitleSize(.14);
    PullHistogram->GetYaxis()->SetTitleOffset(.93);
    PullHistogram->GetYaxis()->SetLabelSize(.12);
  }
  PullHistogram->SetTitle("");
  PullHistogram->SetStats(0);
  return PullHistogram;
}

inline TH1 *getResidual(TH1 *Histogram1, TH1 *Histogram2,
                        TString Name = "res_") {
  if (Histogram1->GetNbinsX() != Histogram2->GetNbinsX() ||
      Histogram1->GetNbinsY() != Histogram2->GetNbinsY() ||
      Histogram1->GetNbinsZ() != Histogram2->GetNbinsZ()) {
    std::cout << "binning doesn't match" << std::endl;
    return 0;
  }
  TH1 *Residuals = (TH1 *)Histogram1->Clone(Name + Histogram1->GetName());
  Residuals->Reset();

  Histogram2->Scale(Histogram1->Integral() / Histogram2->Integral());
  double Limit = 0;
  for (int i = 1; i <= Histogram1->GetNbinsX(); ++i) {
    for (int j = 1; j <= Histogram1->GetNbinsY(); ++j) {
      for (int k = 1; k <= Histogram1->GetNbinsZ(); ++k) {
        unsigned int Bin = Histogram1->GetBin(i, j, k);
        double BinContent1 = Histogram1->GetBinContent(Bin);
        double BinContent2 = Histogram2->GetBinContent(Bin);
        if (BinContent1 == 0 && BinContent2 == 0) {
          double res = -9999;
          Residuals->SetBinContent(Bin, res);
        } else {
          double res = BinContent1 - BinContent2;
          Residuals->SetBinContent(Bin, res);
          if (std::fabs(res) > Limit)
            Limit = std::abs(res);
        }
        Residuals->SetBinError(Bin, 0);
      }
    }
  }
  Limit += 1;
  if (Residuals->GetDimension() == 1)
    Residuals->GetYaxis()->SetRangeUser(
        -Limit, Limit); // symmetric range including all values
  if (Residuals->GetDimension() == 2)
    Residuals->GetZaxis()->SetRangeUser(
        -Limit, Limit); // symmetric range including all values
  return Residuals;
}

inline TPad *drawHist(std::vector<TH1D *> Histograms,
                      std::vector<TString> DrawOption, double Min = 0,
                      double Max = 0) {
  if (Histograms.size() != DrawOption.size())
    throw std::runtime_error("drawPull() | Number of histograms and number "
                             "of draw options does not match!");

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0); // entries only

  TPad *Pad = new TPad();
  if (!Histograms.size())
    return 0;
  for (unsigned int i = 0; i < Histograms.size(); i++)
    Histograms.at(i)->Draw(DrawOption.at(i));
  return Pad;
}

inline TPad *drawPull(std::vector<TH1D *> Histogram,
                      std::vector<TString> DrawOptions, double Min = 0,
                      double Max = 0) {
  if (Histogram.size() != DrawOptions.size())
    throw std::runtime_error("drawPull() | Number of histograms and number "
                             "of draw options does not match!");

  if (!Histogram.size())
    LOG(ERROR) << "drawPull() | No histograms given.";

  Int_t OldOptTitle = gStyle->GetOptTitle();
  Int_t OldOptStat = gStyle->GetOptStat();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0); // entries only

  TPad *Pad = nullptr;
  if (Histogram.size() == 1) {
    Pad = new TPad();
    Histogram.at(0)->Draw(DrawOptions.at(0));
  } else {
    Pad = new TPad("hist", "hist", 0.0, 0.3, 1, 1);
    Pad->Draw();
    Pad->SetMargin(0.1, 0.05, 0.0, 0.05);
    auto PullPad = new TPad("pull", "pull", 0.0, 0.0, 1, 0.3);
    PullPad->Draw();
    PullPad->SetMargin(0.1, 0.05, 0.3, 0.0); // left-right-bottom-top

    if (Histogram.at(0)->GetDimension() != 1 ||
        Histogram.at(1)->GetDimension() != 1) {
      std::cout << "Dimension of histograms larger 1!" << std::endl;
      return Pad;
    }
    if (Histogram.at(0)->GetNbinsX() != Histogram.at(1)->GetNbinsX()) {
      std::cout << "binning doesn't match" << std::endl;
      return Pad;
    }

    Pad->cd();
    Histogram.at(0)->GetYaxis()->SetRangeUser(
        0.00000001, Histogram.at(0)->GetMaximum() * 1.3);
    Histogram.at(0)->GetYaxis()->SetTitleOffset(1.0);
    Histogram.at(0)->GetYaxis()->SetTitleSize(0.05);
    Histogram.at(0)->GetYaxis()->SetLabelSize(0.05);

    Histogram.at(0)->Draw(DrawOptions.at(0));
    Histogram.at(1)->Draw(DrawOptions.at(1));
    for (unsigned int i = 2; i < Histogram.size(); ++i)
      Histogram.at(i)->Draw(DrawOptions.at(i));

    Double_t chi2;
    Double_t Chi2Results[Histogram.at(0)->GetNbinsX()];
    Int_t NDF, IGood;
    Histogram.at(0)->Chi2TestX(Histogram.at(1), chi2, NDF, IGood, "UW",
                               Chi2Results);

    char Chi2Char[60];
    sprintf(Chi2Char, "#chi^{2}_{1D}/ndf = %.2f/%d", chi2, NDF);
    // TLatex* ltx = new TLatex();
    // ltx->DrawLatexNDC(0.2,0.85,chi2Char);//needs to be drawn inside a TPad
    // delete ltx;

    TAxis *XAxis = ((TH1 *)Histogram.at(0))->GetXaxis();
    Int_t NPoints = XAxis->GetNbins();
    Double_t Points[NPoints];
    for (Int_t i = 0; i < NPoints; i++)
      Points[i] = XAxis->GetBinCenter(i + 1);

    TGraph *PullGraph = new TGraph(NPoints, Points, Chi2Results);
    PullGraph->SetName("dataAdaptiveBinned");
    PullPad->cd();

    //	pad_pull->SetGridy();
    PullGraph->Draw("AP"); // draw markers only
    PullGraph->SetTitle("Normalized residuals");
    PullGraph->GetXaxis()->SetTitle(XAxis->GetTitle());
    PullGraph->GetXaxis()->SetTitleOffset(0.8);
    PullGraph->GetXaxis()->SetTitleSize(0.15);
    PullGraph->GetXaxis()->SetLabelSize(0.10);
    PullGraph->GetXaxis()->SetLimits(XAxis->GetBinLowEdge(1),
                                     XAxis->GetBinLowEdge(NPoints + 1));
    PullGraph->GetYaxis()->SetTitle("Dev. [#sigma]");
    PullGraph->GetYaxis()->SetTitleOffset(0.4);
    PullGraph->GetYaxis()->SetNdivisions(504);
    if (Min != Max)
      PullGraph->GetYaxis()->SetRangeUser(Min, Max);
    else
      PullGraph->GetYaxis()->SetRangeUser(-5, 5.5);
    PullGraph->GetYaxis()->SetTitleSize(0.10);
    PullGraph->GetYaxis()->SetLabelSize(0.12);
    PullGraph->GetYaxis()->CenterTitle(1);

    TLine *Line = new TLine(XAxis->GetBinLowEdge(1), 0.0,
                            XAxis->GetBinUpEdge(NPoints + 1), 0.0);
    Line->SetLineStyle(2);
    Line->Draw();
    gPad->Update();

    gStyle->SetOptTitle(OldOptTitle);
    gStyle->SetOptTitle(OldOptStat);
  }
  return Pad;
}

inline TPad *drawPull(TH1D *Histogram1, TH1D *Histogram2, TString DrawOption1,
                      TString DrawOption2, double Min = 0, double Max = 0) {
  std::vector<TH1D *> Histograms;
  Histograms.push_back(Histogram1);
  Histograms.push_back(Histogram2);
  std::vector<TString> DrawOptions;
  DrawOptions.push_back(DrawOption1);
  DrawOptions.push_back(DrawOption2);
  return drawPull(Histograms, DrawOptions, Min, Max);
}

inline TPad *drawResidual(std::vector<TH1D *> Histograms,
                          std::vector<TString> DrawOptions, double Min = 0,
                          double Max = 0) {
  if (Histograms.size() != DrawOptions.size())
    throw std::runtime_error("drawResidual() | Number of histograms and number "
                             "of draw options does not match!");

  Int_t OldOptTitle = gStyle->GetOptTitle();
  Int_t OldOptStat = gStyle->GetOptStat();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0); // entries only

  TPad *Pad = nullptr;
  if (Histograms.size() == 1) {
    Pad = new TPad();
    Histograms.at(0)->Draw(DrawOptions.at(0));
    Histograms.at(0)->GetYaxis()->SetTitleSize(0.06);
  } else if (Histograms.size() == 2) {
    Pad = new TPad("hist", "hist", 0.0, 0.3, 1, 1);
    Pad->Draw();
    Pad->SetMargin(0.1, 0.05, 0.0, 0.05);
    TPad *PullPad = new TPad("pull", "pull", 0.0, 0.0, 1, 0.3);
    PullPad->Draw();
    PullPad->SetMargin(0.1, 0.05, 0.3, 0.0); // left-right-bottom-top

    if (Histograms.at(0)->GetDimension() != 1 ||
        Histograms.at(1)->GetDimension() != 1) {
      std::cout << "Dimension of histograms larger 1!" << std::endl;
      return Pad;
    }
    if (Histograms.at(0)->GetNbinsX() != Histograms.at(1)->GetNbinsX()) {
      std::cout << "binning doesn't match" << std::endl;
      return Pad;
    }

    Pad->cd();
    Histograms.at(0)->GetYaxis()->SetRangeUser(
        0.00000001, Histograms.at(0)->GetMaximum() * 1.1);
    Histograms.at(0)->GetYaxis()->SetTitleOffset(0.83);

    Histograms.at(0)->Draw(DrawOptions.at(0));
    Histograms.at(1)->Draw(DrawOptions.at(1));
    Histograms.at(0)->GetYaxis()->SetTitleSize(0.06);
    for (unsigned int i = 2; i < Histograms.size(); ++i) {
      Histograms.at(i)->Draw(DrawOptions.at(i));
    }

    Double_t Chi2;
    Double_t Chi2Result[Histograms.at(0)->GetNbinsX()];
    Int_t NDF, IGood;
    Histograms.at(0)->Chi2TestX(Histograms.at(1), Chi2, NDF, IGood, "UW",
                                Chi2Result);

    char Chi2Char[60];
    sprintf(Chi2Char, "#chi^{2}_{1D}/ndf = %.2f/%d", Chi2, NDF);

    TAxis *XAxis = ((TH1 *)Histograms.at(0))->GetXaxis();
    Int_t NPoints = XAxis->GetNbins();

    TH1D *Residuals = (TH1D *)Histograms.at(0)->Clone("resHist");
    Residuals->Add(Histograms.at(1), -1);
    PullPad->cd();
    Residuals->Draw("E"); // draw markers only
    Residuals->SetTitle("Residuals");
    Residuals->GetXaxis()->SetTitle(XAxis->GetTitle());
    Residuals->GetXaxis()->SetTitleOffset(0.8);
    Residuals->GetXaxis()->SetTitleSize(0.15);
    Residuals->GetXaxis()->SetLabelSize(0.10);
    Residuals->GetXaxis()->SetLimits(XAxis->GetBinLowEdge(1),
                                     XAxis->GetBinLowEdge(NPoints + 1));
    Residuals->GetYaxis()->SetTitle("Deviation");
    Residuals->GetYaxis()->SetTitleOffset(0.4);
    Residuals->GetYaxis()->SetNdivisions(504);
    if (Min != Max)
      Residuals->GetYaxis()->SetRangeUser(Min, Max);
    Residuals->GetYaxis()->SetTitleSize(0.12);
    Residuals->GetYaxis()->SetLabelSize(0.12);
    Residuals->GetYaxis()->CenterTitle(1);
    Residuals->SetMarkerStyle(20);

    TLine *Line = new TLine(XAxis->GetBinLowEdge(1), 0.0,
                            XAxis->GetBinUpEdge(NPoints + 1), 0.0);
    Line->SetLineStyle(2);
    Line->Draw();
    gPad->Update();

    gStyle->SetOptTitle(OldOptTitle);
    gStyle->SetOptTitle(OldOptStat);
  }
  return Pad;
}

inline TPad *drawResidual(TH1D *Histogram1, TH1D *Histogram2,
                          TString DrawOption1, TString DrawOption2,
                          double Min = 0, double Max = 0) {
  std::vector<TH1D *> Histograms;
  Histograms.push_back(Histogram1);
  Histograms.push_back(Histogram2);
  std::vector<TString> DrawOptions;
  DrawOptions.push_back(DrawOption1);
  DrawOptions.push_back(DrawOption2);
  return drawResidual(Histograms, DrawOptions, Min, Max);
}

inline void getTH2PolyChi2(TH2Poly *Histogram1, TH2Poly *Histogram2,
                           double &Chi2, int &NDF, int &IGood) {
  if (Histogram1->GetBins()->GetEntries() !=
      Histogram2->GetBins()->GetEntries()) {
    std::cout << "binning doesn't match" << std::endl;
    return;
  }
  Chi2 = 0;
  NDF = 0;
  for (int Bin = 1; Bin <= Histogram1->GetBins()->GetEntries(); Bin++) {
    double BinContent1 = Histogram1->GetBinContent(Bin);
    double BinContent2 = Histogram2->GetBinContent(Bin);
    double BinError1 = Histogram1->GetBinError(Bin);
    double BinError2 = Histogram2->GetBinError(Bin);
    if (BinContent1 > 0 || BinContent2 > 0) {
      Chi2 += (BinContent1 - BinContent2) * (BinContent1 - BinContent2) /
              (BinError1 * BinError1 + BinError2 * BinError2);
      NDF++;
    }
  }
  return;
}

inline TH2Poly *getTH2PolyPull(TH2Poly *Histogram1, TH2Poly *Histogram2,
                               TString Name) {
  if (Histogram1->GetBins()->GetEntries() !=
      Histogram2->GetBins()->GetEntries()) {
    std::cout << "binning doesn't match" << std::endl;
    return 0;
  }
  TH2Poly *Residuals =
      (TH2Poly *)Histogram1->Clone(Name + Histogram1->GetName());
  Residuals->Reset("");

  Histogram2->Scale(Histogram1->Integral() / Histogram2->Integral());
  double Limit = 0;
  double Integral = 0;
  int NBins = 0;
  for (int bin = 1; bin <= Histogram1->GetBins()->GetEntries(); bin++) {
    double BinContent1 = Histogram1->GetBinContent(bin);
    double BinContent2 = Histogram2->GetBinContent(bin);
    double BinError1 = Histogram1->GetBinError(bin);
    double BinError2 = Histogram2->GetBinError(bin);
    if (BinContent1 > 0 && BinContent2 > 0) {
      double Pull = (BinContent1 - BinContent2) /
                    sqrt(BinError1 * BinError1 + BinError2 * BinError2);
      Integral += Pull;
      NBins++;
      Residuals->SetBinContent(bin, Pull);
      if (std::fabs(Pull) > Limit)
        Limit = std::abs(Pull);
    } else
      Residuals->SetBinContent(bin, -999); // empty bins are set to error value
    Residuals->SetBinError(bin, 0);
  }
  Limit = 4;
  Limit += 1;
  Residuals->GetZaxis()->SetRangeUser(
      -Limit, Limit); // symmetric range including all values
  Residuals->GetZaxis()->SetTitle("deviation [#sigma]");
  Residuals->GetZaxis()->SetTitleOffset(0.5);
  std::cout << "Creating pull histogram for TH2Poly! integral=" << Integral
            << " nBins=" << NBins << std::endl;

  return Residuals;
}

inline TH2Poly *adaptiveBinning(UInt_t DataSize, UInt_t DataDimension,
                                Double_t *Data, UInt_t NBins = 100) {
  UInt_t Size =
      UInt_t(DataSize / NBins) * NBins; // size should be a multiple of nBins
  if (Size == 0) {                      // return empty histogram
    TH2Poly *Histogram = new TH2Poly();
    Histogram->AddBin(0, 1, 0, 1);
    return Histogram;
  }

  TKDTreeBinning KDBins(Size, DataDimension, Data, NBins);

  UInt_t NKDBins = KDBins.GetNBins();
  UInt_t Dimension = KDBins.GetDim();

  const Double_t *binsMinEdges = KDBins.GetBinsMinEdges();
  const Double_t *binsMaxEdges = KDBins.GetBinsMaxEdges();

  TH2Poly *Histogram = new TH2Poly("h2PolyBinTest", "", KDBins.GetDataMin(0),
                                   KDBins.GetDataMax(0), KDBins.GetDataMin(1),
                                   KDBins.GetDataMax(1));

  for (UInt_t i = 0; i < NKDBins; ++i) {
    UInt_t EdgeDim = i * Dimension;
    Histogram->AddBin(binsMinEdges[EdgeDim], binsMinEdges[EdgeDim + 1],
                      binsMaxEdges[EdgeDim], binsMaxEdges[EdgeDim + 1]);
  }

  return Histogram;
}
#endif
