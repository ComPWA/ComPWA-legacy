//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <memory>
#include "PolyFit.hpp"

#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TF1.h"
#include "TGraphErrors.h"

//#include "ErrLogger/ErrLogger.hh"

std::shared_ptr<ControlParameter> PolyFit::createInstance(double p0, double p1, double p2, double p3, double sigma){
  if(!instance_)
    instance_ = std::shared_ptr<ControlParameter>(new PolyFit(p0, p1, p2, p3, sigma));

  return instance_;
}

PolyFit::PolyFit(double p0, double p1, double p2, double p3, double sigma) :
  _theTFile(new TFile("test/myFit.root","RECREATE")),
  _sigma(sigma)
{

  // Display parameters for test distribution 

  std::cout << std::endl;
  /*Info <<"Set p0 as "<< p0 << endmsg;
  Info <<"Set p1 as "<< p1 << endmsg;
  Info <<"Set p2 as "<< p2 << endmsg;
  Info <<"Set p3 as "<< p3 << endmsg;
  Info <<"Set sigma as "<< sigma << endmsg;*/
  std::cout <<"Set p0 as "<< p0 << std::endl;
  std::cout <<"Set p1 as "<< p1 << std::endl;
  std::cout <<"Set p2 as "<< p2 << std::endl;
  std::cout <<"Set p3 as "<< p3 << std::endl;
  std::cout <<"Set sigma as "<< sigma << std::endl;

  // Generate test distribution and smear them with a gaussian
  TRandom randomNumber;
  for(int i=0; i<1000; i++) {
    double tmpXvalue=static_cast<double>((rand() % 10000 + 1)) / 100;
    double tmpYvalue=randomNumber.Gaus(p0 + p1 * tmpXvalue + p2 * tmpXvalue * tmpXvalue + p3 * tmpXvalue * tmpXvalue * tmpXvalue, _sigma);
    _xValue.push_back(tmpXvalue);
    _yValue.push_back(tmpYvalue);
  }

}

double PolyFit::controlParameter(ParameterList& minPar){
 
  // Calculate chi^2 for current set of fit parameters 
  double result=0.;
  for (unsigned int i=0; i<_xValue.size(); i++){
    double yValFit=minPar.GetParameterValue(0)+minPar.GetParameterValue(1)*_xValue[i]+minPar.GetParameterValue(2)*_xValue[i]*_xValue[i]+minPar.GetParameterValue(3)*_xValue[i]*_xValue[i]*_xValue[i];
    double yValExp=_yValue[i];
    double tmpChi=((yValExp-yValFit)*(yValExp-yValFit))/(_sigma*_sigma);
    result+=tmpChi;
  }
  return result;
}


void PolyFit::drawGraph(double a, double b, double c, double d){

  // Create arrays (which will be needed for feeding the TGraph)
  const unsigned int dataEntries=_xValue.size();
  double x[dataEntries];
  double y[dataEntries]; 
  double xerr[dataEntries];
  double yerr[dataEntries]; 

  // Convert vectors to arrays
  for(unsigned int i=0; i<dataEntries; i++) {
    x[i]=_xValue[i];
    y[i]=_yValue[i];
    xerr[i] = 0;
    yerr[i] = _sigma;
  }

  // Create ROOT file                                                                                                             
  //_theTFile = new TFile("myFit.root","RECREATE");

  // Create TGraph and write to ROOT file 
  TCanvas* c1=new TCanvas("c1","c1",1200,800);
  c1->cd();
  TGraphErrors* linDist = new TGraphErrors(100, x, y, xerr, yerr);
  linDist->SetMarkerStyle(6);
  linDist->Draw("AP");
  linDist->Write("myGraph");

  // Draw fit with calculated parameters
  TF1* fit1 = new TF1("fit1","pol3",0.,100.);
  fit1->SetParameters(a, b, c, d);
  fit1->SetLineColor(46);
  fit1->SetLineWidth(3);
  fit1->Draw();
  fit1->Write();

}


PolyFit::~PolyFit()
{
  //cout << "writa polyy-fit data" << endl;
  _theTFile->Write();
  _theTFile->Close();
  //cout << "done" << endl << endl;
  //delete instance_;
}

