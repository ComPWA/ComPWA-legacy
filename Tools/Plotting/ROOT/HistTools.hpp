// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Collection of hopefully usefule routines to plot histograms,
/// e.g. creating (normalized) residual plots
///

#ifndef TOOLS_H
#define TOOLS_H

#include <TGraph.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TPad.h>
#include <TMath.h>
#include <TKDTreeBinning.h>


///
///
///
inline TH1 *getPull(TH1 *hist1, TH1 *hist2, TString name = "pull_") {
  if (hist1->GetNbinsX() != hist2->GetNbinsX() ||
      hist1->GetNbinsY() != hist2->GetNbinsY() ||
      hist1->GetNbinsZ() != hist2->GetNbinsZ()) {
    std::cout << "binning doesnt match" << std::endl;
    return 0;
  }
  TH1 *pullHist = (TH1 *)hist1->Clone(name + hist1->GetName());
  pullHist->Reset();
  double limitPull = 0;
  for (int i = 1; i <= hist1->GetNbinsX(); ++i) {
    for (int j = 1; j <= hist1->GetNbinsY(); ++j) {
      for (int k = 1; k <= hist1->GetNbinsZ(); ++k) {
        unsigned int bin = hist1->GetBin(i, j, k);
        double c1 = hist1->GetBinContent(bin);
        double c2 = hist2->GetBinContent(bin);
        double c1Err = hist1->GetBinError(bin);
        double c2Err = hist2->GetBinError(bin);
        double pull;
        if (c1 > 0 && c2 > 0) {
          pull = (c1 - c2) / sqrt(c1Err * c1Err + c2Err * c2Err);
          pullHist->SetBinContent(bin, pull);
          if (std::fabs(pull) > limitPull)
            limitPull = std::abs(pull);
        } else
          pullHist->SetBinContent(bin,
                                  -999); // empty bins are set to error value
        pullHist->SetBinError(bin, .0001);
      }
    }
  }
  limitPull = 4; // set fix limits
  limitPull += 1;
  // X
  pullHist->GetXaxis()->SetTitleSize(.14);
  pullHist->GetXaxis()->SetTitleOffset(.93);
  pullHist->GetXaxis()->SetLabelSize(.12);
  if (pullHist->GetDimension() == 1) {
    pullHist->GetYaxis()->SetRangeUser(-limitPull, limitPull);
    pullHist->GetYaxis()->SetTitle("deviation [#sigma]");
    pullHist->GetYaxis()->SetTitleOffset(0.36);
    pullHist->GetYaxis()->SetTitleSize(0.12);
    pullHist->GetYaxis()->SetLabelSize(0.12);
    pullHist->GetYaxis()->SetNdivisions(504);
    pullHist->GetYaxis()->CenterTitle();
  }
  if (pullHist->GetDimension() == 2) { // symmetric range including all values
    pullHist->GetZaxis()->SetRangeUser(-limitPull, limitPull);
    pullHist->GetZaxis()->SetTitle("deviation [#sigma]");
    pullHist->GetYaxis()->SetTitleSize(.14);
    pullHist->GetYaxis()->SetTitleOffset(.93);
    pullHist->GetYaxis()->SetLabelSize(.12);
  }
  pullHist->SetTitle("");
  pullHist->SetStats(0);
  return pullHist;
}

///
///
///
inline TH1 *getResidual(TH1 *hist1, TH1 *hist2, TString name = "res_") {
  if (hist1->GetNbinsX() != hist2->GetNbinsX() ||
      hist1->GetNbinsY() != hist2->GetNbinsY() ||
      hist1->GetNbinsZ() != hist2->GetNbinsZ()) {
    std::cout << "binning doesnt match" << std::endl;
    return 0;
  }
  TH1 *resHist = (TH1 *)hist1->Clone(name + hist1->GetName());
  resHist->Reset();

  hist2->Scale(hist1->Integral() / hist2->Integral());
  double limit = 0;
  for (int i = 1; i <= hist1->GetNbinsX(); ++i) {
    for (int j = 1; j <= hist1->GetNbinsY(); ++j) {
      for (int k = 1; k <= hist1->GetNbinsZ(); ++k) {
        unsigned int bin = hist1->GetBin(i, j, k);
        double c1 = hist1->GetBinContent(bin);
        double c2 = hist2->GetBinContent(bin);
        if (c1 == 0 && c2 == 0) {
          double res = -9999;
          resHist->SetBinContent(bin, res);
        } else {
          double res = c1 - c2;
          resHist->SetBinContent(bin, res);
          if (std::fabs(res) > limit)
            limit = std::abs(res);
        }
        resHist->SetBinError(bin, 0);
      }
    }
  }
  limit += 1;
  if (resHist->GetDimension() == 1)
    resHist->GetYaxis()->SetRangeUser(
        -limit, limit); // symmetric range including all values
  if (resHist->GetDimension() == 2)
    resHist->GetZaxis()->SetRangeUser(
        -limit, limit); // symmetric range including all values
  return resHist;
}

///
///
///
inline TPad *drawHist(std::vector<TH1D *> hist, std::vector<TString> drawOption,
                      double min = 0, double max = 0) {
  if (hist.size() != drawOption.size())
    throw std::runtime_error("drawPull() | Number of histograms and number "
                             "of draw options does not match!");

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0); // entries only

  TPad *pad = new TPad();
  if (!hist.size())
    return 0;
  for (unsigned int i = 0; i < hist.size(); i++)
    hist.at(i)->Draw(drawOption.at(i));
  return pad;
}

///
///
///
inline TPad *drawPull(std::vector<TH1D *> hist, std::vector<TString> drawOption,
                      double min = 0, double max = 0) {
  if (hist.size() != drawOption.size())
    throw std::runtime_error("drawPull() | Number of histograms and number "
                             "of draw options does not match!");

  if (!hist.size())
    LOG(ERROR) << "drawPull() | No histograms given.";
    
  Int_t optTitle = gStyle->GetOptTitle();
  Int_t optStat = gStyle->GetOptStat();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0); // entries only

  TPad *pad = NULL;
  if (hist.size() == 1) {
    pad = new TPad();
    hist.at(0)->Draw(drawOption.at(0));
  } else {
    pad = new TPad("hist", "hist", 0.0, 0.3, 1, 1);
    pad->Draw();
    pad->SetMargin(0.1, 0.05, 0.0, 0.05);
    TPad *pad_pull = new TPad("pull", "pull", 0.0, 0.0, 1, 0.3);
    pad_pull->Draw();
    pad_pull->SetMargin(0.1, 0.05, 0.3, 0.0); // left-right-bottom-top

    if (hist.at(0)->GetDimension() != 1 || hist.at(1)->GetDimension() != 1) {
      std::cout << "Dimension of histograms larger 1!" << std::endl;
      return pad;
    }
    if (hist.at(0)->GetNbinsX() != hist.at(1)->GetNbinsX()) {
      std::cout << "binning doesnt match" << std::endl;
      return pad;
    }

    pad->cd();
    hist.at(0)->GetYaxis()->SetRangeUser(0.00000001,
                                         hist.at(0)->GetMaximum() * 1.3);
    hist.at(0)->GetYaxis()->SetTitleOffset(0.83);

    hist.at(0)->Draw(drawOption.at(0));
    hist.at(1)->Draw(drawOption.at(1));
    for (unsigned int i = 2; i < hist.size(); ++i)
      hist.at(i)->Draw(drawOption.at(i));

    Double_t chi2;
    Double_t res[hist.at(0)->GetNbinsX()];
    Int_t ndf, igood;
    hist.at(0)->Chi2TestX(hist.at(1), chi2, ndf, igood, "UW", res);

    char chi2Char[60];
    sprintf(chi2Char, "#chi^{2}_{1D}/ndf = %.2f/%d", chi2, ndf);
    // TLatex* ltx = new TLatex();
    // ltx->DrawLatexNDC(0.2,0.85,chi2Char);//needs to be drawn inside a TPad
    // delete ltx;

    TAxis *xaxis = ((TH1 *)hist.at(0))->GetXaxis();
    Int_t fNpoints = xaxis->GetNbins();
    Double_t fX[fNpoints];
    for (Int_t i = 0; i < fNpoints; i++)
      fX[i] = xaxis->GetBinCenter(i + 1);

    TGraph *pullGr = new TGraph(fNpoints, fX, res);
    pullGr->SetName("dataAdaptiveBinned");
    pad_pull->cd();

    //	pad_pull->SetGridy();
    pullGr->Draw("AP"); // draw markers only
    pullGr->SetTitle("Normalized residuals");
    pullGr->GetXaxis()->SetTitle(xaxis->GetTitle());
    pullGr->GetXaxis()->SetTitleOffset(0.8);
    pullGr->GetXaxis()->SetTitleSize(0.15);
    pullGr->GetXaxis()->SetLabelSize(0.10);
    pullGr->GetXaxis()->SetLimits(xaxis->GetBinLowEdge(1),
                                  xaxis->GetBinLowEdge(fNpoints + 1));
    pullGr->GetYaxis()->SetTitle("Deviation [#sigma]");
    pullGr->GetYaxis()->SetTitleOffset(0.4);
    pullGr->GetYaxis()->SetNdivisions(504);
    if (min != max)
      pullGr->GetYaxis()->SetRangeUser(min, max);
    else
      pullGr->GetYaxis()->SetRangeUser(-5, 5.5);
    pullGr->GetYaxis()->SetTitleSize(0.12);
    pullGr->GetYaxis()->SetLabelSize(0.12);
    pullGr->GetYaxis()->CenterTitle(1);

    TLine *line = new TLine(xaxis->GetBinLowEdge(1), 0.0,
                            xaxis->GetBinUpEdge(fNpoints + 1), 0.0);
    line->SetLineStyle(2);
    line->Draw();
    gPad->Update();

    gStyle->SetOptTitle(optTitle);
    gStyle->SetOptTitle(optStat);
  }
  return pad;
}

///
///
///
inline TPad *drawPull(TH1D *hist1, TH1D *hist2, TString drawOption1,
                      TString drawOption2, double min = 0, double max = 0) {
  std::vector<TH1D *> vHist;
  vHist.push_back(hist1);
  vHist.push_back(hist2);
  std::vector<TString> vOpt;
  vOpt.push_back(drawOption1);
  vOpt.push_back(drawOption2);
  return drawPull(vHist, vOpt, min, max);
}

///
///
///
inline TPad *drawResidual(std::vector<TH1D *> hist,
                          std::vector<TString> drawOption, double min = 0,
                          double max = 0) {
  if (hist.size() != drawOption.size())
    throw std::runtime_error("drawResidual() | Number of histograms and number "
                             "of draw options does not match!");

  Int_t optTitle = gStyle->GetOptTitle();
  Int_t optStat = gStyle->GetOptStat();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0); // entries only

  TPad *pad = NULL;
  if (hist.size() == 1) {
    pad = new TPad();
    hist.at(0)->Draw(drawOption.at(0));
  } else if (hist.size() == 2) {
    pad = new TPad("hist", "hist", 0.0, 0.3, 1, 1);
    pad->Draw();
    pad->SetMargin(0.1, 0.05, 0.0, 0.05);
    TPad *pad_pull = new TPad("pull", "pull", 0.0, 0.0, 1, 0.3);
    pad_pull->Draw();
    pad_pull->SetMargin(0.1, 0.05, 0.3, 0.0); // left-right-bottom-top

    if (hist.at(0)->GetDimension() != 1 || hist.at(1)->GetDimension() != 1) {
      std::cout << "Dimension of histograms larger 1!" << std::endl;
      return pad;
    }
    if (hist.at(0)->GetNbinsX() != hist.at(1)->GetNbinsX()) {
      std::cout << "binning doesnt match" << std::endl;
      return pad;
    }

    pad->cd();
    hist.at(0)->GetYaxis()->SetRangeUser(0.00000001,
                                         hist.at(0)->GetMaximum() * 1.1);
    hist.at(0)->GetYaxis()->SetTitleOffset(0.83);

    hist.at(0)->Draw(drawOption.at(0));
    hist.at(1)->Draw(drawOption.at(1));
    for (unsigned int i = 2; i < hist.size(); ++i)
      hist.at(i)->Draw(drawOption.at(i));

    Double_t chi2;
    Double_t res[hist.at(0)->GetNbinsX()];
    Int_t ndf, igood;
    hist.at(0)->Chi2TestX(hist.at(1), chi2, ndf, igood, "UW", res);

    char chi2Char[60];
    sprintf(chi2Char, "#chi^{2}_{1D}/ndf = %.2f/%d", chi2, ndf);
    // TLatex* ltx = new TLatex();
    // ltx->DrawLatexNDC(0.2,0.85,chi2Char);//needs to be drawn inside a TPad
    // delete ltx;

    TAxis *xaxis = ((TH1 *)hist.at(0))->GetXaxis();
    Int_t fNpoints = xaxis->GetNbins();
    //		Double_t fX[fNpoints];
    //		for (Int_t i = 0; i < fNpoints; i++)
    //			fX[i] = xaxis->GetBinCenter(i + 1);

    TH1D *resHist = (TH1D *)hist.at(0)->Clone("resHist");
    resHist->Add(hist.at(1), -1);
    pad_pull->cd();
    resHist->Draw("E"); // draw markers only
    resHist->SetTitle("Residuals");
    resHist->GetXaxis()->SetTitle(xaxis->GetTitle());
    resHist->GetXaxis()->SetTitleOffset(0.8);
    resHist->GetXaxis()->SetTitleSize(0.15);
    resHist->GetXaxis()->SetLabelSize(0.10);
    resHist->GetXaxis()->SetLimits(xaxis->GetBinLowEdge(1),
                                   xaxis->GetBinLowEdge(fNpoints + 1));
    resHist->GetYaxis()->SetTitle("Deviation");
    resHist->GetYaxis()->SetTitleOffset(0.4);
    resHist->GetYaxis()->SetNdivisions(504);
    if (min != max)
      resHist->GetYaxis()->SetRangeUser(min, max);
    resHist->GetYaxis()->SetTitleSize(0.12);
    resHist->GetYaxis()->SetLabelSize(0.12);
    resHist->GetYaxis()->CenterTitle(1);

    TLine *line = new TLine(xaxis->GetBinLowEdge(1), 0.0,
                            xaxis->GetBinUpEdge(fNpoints + 1), 0.0);
    line->SetLineStyle(2);
    line->Draw();
    gPad->Update();

    gStyle->SetOptTitle(optTitle);
    gStyle->SetOptTitle(optStat);
  }
  return pad;
}

///
///
///
inline TPad *drawResidual(TH1D *hist1, TH1D *hist2, TString drawOption1,
                          TString drawOption2, double min = 0, double max = 0) {
  std::vector<TH1D *> vHist;
  vHist.push_back(hist1);
  vHist.push_back(hist2);
  std::vector<TString> vOpt;
  vOpt.push_back(drawOption1);
  vOpt.push_back(drawOption2);
  return drawResidual(vHist, vOpt, min, max);
}

///
///
///
inline void getTH2PolyChi2(TH2Poly *hist1, TH2Poly *hist2, double &chi2,
                           int &ndf, int &igood) {
  if (hist1->GetBins()->GetEntries() != hist2->GetBins()->GetEntries()) {
    std::cout << "binning doesnt match" << std::endl;
    return;
  }
  //	std::cout<<hist1->GetEntries()<<std::endl;
  //	std::cout<<hist1->GetIntegral()<<std::endl;
  //	std::cout<<hist1->GetSumOfWeights()<<std::endl;
  chi2 = 0;
  ndf = 0;
  for (int bin = 1; bin <= hist1->GetBins()->GetEntries(); bin++) {
    double c1 = hist1->GetBinContent(bin);
    double c2 = hist2->GetBinContent(bin);
    double c1Err = hist1->GetBinError(bin);
    double c2Err = hist2->GetBinError(bin);
    if (c1 > 0 || c2 > 0) {
      chi2 += (c1 - c2) * (c1 - c2) / (c1Err * c1Err + c2Err * c2Err);
      //			chi2+=(c1-c2)*(c1-c2)/(c1Err*c1Err);
      //			std::cout<<chi2<<std::endl;
      ndf++;
    }
  }
  //	std::cout<<"Calculating TH2Poly chi2="<<chi2<<" ndf="<<ndf<<std::endl;
  return;
}

///
///
///
inline TH2Poly *getTH2PolyPull(TH2Poly *hist1, TH2Poly *hist2, TString name) {
  if (hist1->GetBins()->GetEntries() != hist2->GetBins()->GetEntries()) {
    std::cout << "binning doesnt match" << std::endl;
    return 0;
  }
  TH2Poly *resHist = (TH2Poly *)hist1->Clone(name + hist1->GetName());
  resHist->Reset("");

  hist2->Scale(hist1->Integral() / hist2->Integral());
  double limit = 0;
  double integral = 0;
  int nBins = 0;
  for (int bin = 1; bin <= hist1->GetBins()->GetEntries(); bin++) {
    double c1 = hist1->GetBinContent(bin);
    double c2 = hist2->GetBinContent(bin);
    double c1Err = hist1->GetBinError(bin);
    double c2Err = hist2->GetBinError(bin);
    if (c1 > 0 && c2 > 0) {
      double pull = (c1 - c2) / sqrt(c1Err * c1Err + c2Err * c2Err);
      integral += pull;
      nBins++;
      resHist->SetBinContent(bin, pull);
      if (std::fabs(pull) > limit)
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
  std::cout << "Creating pull histogram for TH2Poly! integral=" << integral
            << " nBins=" << nBins << std::endl;

  return resHist;
}

///
///
///
inline TH2Poly *adaptiveBinning(UInt_t dataSize, UInt_t dataDim, Double_t *data,
                                UInt_t nBins = 100) {
  UInt_t size =
      UInt_t(dataSize / nBins) * nBins; // size should be a multiple of nBins
  if (size == 0) {                      // return empty histogram
    TH2Poly *h2p = new TH2Poly();
    h2p->AddBin(0, 1, 0, 1);
    return h2p;
  }

  TKDTreeBinning kdBins(size, dataDim, data, nBins);

  UInt_t nbins = kdBins.GetNBins();
  UInt_t dim = kdBins.GetDim();

  const Double_t *binsMinEdges = kdBins.GetBinsMinEdges();
  const Double_t *binsMaxEdges = kdBins.GetBinsMaxEdges();

  TH2Poly *h2pol = new TH2Poly("h2PolyBinTest", "", kdBins.GetDataMin(0),
                               kdBins.GetDataMax(0), kdBins.GetDataMin(1),
                               kdBins.GetDataMax(1));

  for (UInt_t i = 0; i < nbins; ++i) {
    UInt_t edgeDim = i * dim;
    h2pol->AddBin(binsMinEdges[edgeDim], binsMinEdges[edgeDim + 1],
                  binsMaxEdges[edgeDim], binsMaxEdges[edgeDim + 1]);
  }

  return h2pol;
}
#endif
