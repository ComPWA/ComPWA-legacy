#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TLatex.h"
#include "TString.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TPad.h"
#include "TFile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <vector>
#include <cstdlib>

using namespace std;

//constants
  const Double_t M = 3.096916; // GeV/c² (J/psi+)
  const Double_t Br = 0.000093; // GeV/c² (width)
  const Double_t m1 = 0.; // GeV/c² (gamma)
  const Double_t m2 = 0.139570; // GeV/c² (pi)
  const Double_t m3 = 0.139570; // GeV/c² (pi)
  //const Double_t c = 299792458.; // m/s
  const Double_t PI = 3.14159; // m/s

  const Double_t m23_sq_min = (m2+m3)*(m2+m3);
  const Double_t m23_sq_max = (M-m1)*(M-m1);
  const Double_t m13_sq_min = (m1+m3)*(m1+m3);
  const Double_t m13_sq_max = (M-m2)*(M-m2);
  const Double_t m12_sq_min = (m1+m2)*(m1+m2);
  const Double_t m12_sq_max = (M-m3)*(M-m3);

  const Double_t m23_min = (m2+m3);
  const Double_t m23_max = (M-m1);
  const Double_t m13_min = (m1+m3);
  const Double_t m13_max = (M-m2);
  const Double_t m12_min = (m1+m2);
  const Double_t m12_max = (M-m3);

double lambda(double x, double y, double z){
  return x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z;
}

double m13_sq_max_constr(Double_t *x, Double_t *par){
  double m23_sq=x[0];
  return m1*m1+m3*m3+0.5/m23_sq*((M*M-m23_sq-m1*m1)*(m23_sq-m2*m2+m3*m3)+sqrt(lambda(m23_sq,M*M,m1*m1))*sqrt(lambda(m23_sq,m2*m2,m3*m3)));
}

double m13_sq_min_constr(Double_t *x, Double_t *par){
  double m23_sq=x[0];
  return m1*m1+m3*m3+0.5/m23_sq*((M*M-m23_sq-m1*m1)*(m23_sq-m2*m2+m3*m3)-sqrt(lambda(m23_sq,M*M,m1*m1))*sqrt(lambda(m23_sq,m2*m2,m3*m3)));
}

double m13_sq_max_constr_varMP(Double_t *x, Double_t *par){
  double m23_sq=x[0];
  double vM=M+3*Br;
  return m1*m1+m3*m3+0.5/m23_sq*((vM*vM-m23_sq-m1*m1)*(m23_sq-m2*m2+m3*m3)+sqrt(lambda(m23_sq,vM*vM,m1*m1))*sqrt(lambda(m23_sq,m2*m2,m3*m3)));
}

double m13_sq_min_constr_varMP(Double_t *x, Double_t *par){
  double m23_sq=x[0];
  double vM=M+3*Br;
  return m1*m1+m3*m3+0.5/m23_sq*((vM*vM-m23_sq-m1*m1)*(m23_sq-m2*m2+m3*m3)-sqrt(lambda(m23_sq,vM*vM,m1*m1))*sqrt(lambda(m23_sq,m2*m2,m3*m3)));
}

double m13_sq_max_constr_varMM(Double_t *x, Double_t *par){
  double m23_sq=x[0];
  double vM=M-3*Br;
  return m1*m1+m3*m3+0.5/m23_sq*((vM*vM-m23_sq-m1*m1)*(m23_sq-m2*m2+m3*m3)+sqrt(lambda(m23_sq,vM*vM,m1*m1))*sqrt(lambda(m23_sq,m2*m2,m3*m3)));
}

double m13_sq_min_constr_varMM(Double_t *x, Double_t *par){
  double m23_sq=x[0];
  double vM=M-3*Br;
  return m1*m1+m3*m3+0.5/m23_sq*((vM*vM-m23_sq-m1*m1)*(m23_sq-m2*m2+m3*m3)-sqrt(lambda(m23_sq,vM*vM,m1*m1))*sqrt(lambda(m23_sq,m2*m2,m3*m3)));
}

double Get_q(double m){
	//const double m0 = 0.493;
	//const double m1 = 0.139;
	//const double m2 = 0.139;
	const double m02 = m1*m1;
	const double m12 = m2*m2;
	const double m22 = m3*m3;
	const double m04 = m02*m02;
	const double m14 = m12*m12;
	const double m24 = m22*m22;
	const double factor = m04+m14+m24-2.*(m02*m12+m12*m22+m22*m02)/4.;
	return sqrt(fabs(factor/m*m));
}

double qValue(double x, double ma=m2, double mb=m3){
	double mapb = ma + mb;
	double mamb = ma - mb;
	double xSq =x*x;
	double t1 = xSq - mapb*mapb;
	double t2 = xSq - mamb*mamb;

	if( t1 < 0 ) {
		//std::cout<<"AmpKinematics: Trying to calculate break-up momentum below threshold!"<<std::endl;
		return 1; //below threshold
	}
	double result=sqrt( t1 * t2 ) / (2. * x );
	return result;

}

double Get_phasespace_factor(double m, double m0){
	double q0 = Get_q(m0);
	double q  = Get_q(m);
	return m0*q/(m*q0);
}

void JakePlotter(TString fileName="JPSIDATA.ACC.root"){

	gROOT->SetStyle("Plain");
	gStyle->SetTitleFont(10*13+2,"xyz");
	gStyle->SetTitleSize(0.06, "xyz");
	gStyle->SetTitleOffset(1.3,"y");
	gStyle->SetTitleOffset(1.3,"z");
	gStyle->SetLabelFont(10*13+2,"xyz");
	gStyle->SetLabelSize(0.06,"xyz");
	gStyle->SetLabelOffset(0.009,"xyz");
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetPadTopMargin(0.16);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadRightMargin(0.16);
	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(0);
	gROOT->ForceStyle();
	gStyle->SetFrameFillColor(0);
   	gStyle->SetFrameFillStyle(0);
   	//TGaxis::SetMaxDigits(3);
	
    double fe[3];
    double fpx[3];
    double fpy[3];
    double fpz[3];
    double feventWeight;


   TFile* fFile = new TFile(fileName);
   TTree* fTree = (TTree*) fFile->Get("kin");

   fTree->SetBranchAddress("e",fe);
   fTree->SetBranchAddress("px",fpx);
   fTree->SetBranchAddress("py",fpy);
   fTree->SetBranchAddress("pz",fpz);
   fTree->SetBranchAddress("weight",&feventWeight);

   unsigned int nParts = 3;
   TH1D *wP[3], *w23 = new TH1D("w23","CosTheta of final state particle 2 & 3",200,-1.,1.);
   w23->GetXaxis()->SetTitle("CosTheta");
   w23->GetYaxis()->SetTitle("Entries");
   for(unsigned int i=0; i<nParts; i++){
     TString name="CosTheta of final state particle ";
     name+=i;
     wP[i] = new TH1D(name,name,200,-1.,1.);
     wP[i]->GetXaxis()->SetTitle("CosTheta");
     wP[i]->GetYaxis()->SetTitle("Entries");
   }

   for(unsigned int evt=0; evt<fTree->GetEntries(); evt++){
	  fTree->GetEntry(evt);

	  vector<TLorentzVector> fPart;

	  for(unsigned int part=0; part<nParts; part++){
		 fPart.push_back(TLorentzVector(fpx[part], fpy[part], fpz[part], fe[part]));
                 wP[part]->Fill(fPart[part].CosTheta());
	  }//particle loop
          w23->Fill(fPart[1].CosTheta());
          w23->Fill(fPart[2].CosTheta());

          if(evt%10000==0){
            cout << "Event " << evt << endl;
            cout << "Particle 0: " << fpx[0] << " " << fpy[0] << " " << fpz[0] << " " << fe[0] << endl;
            cout << "Particle 1: " << fpx[1] << " " << fpy[1] << " " << fpz[1] << " " << fe[1] << endl;
            cout << "Particle 2: " << fpx[2] << " " << fpy[2] << " " << fpz[2] << " " << fe[2] << endl;
            cout << "TLVec 0: " << fPart[0].X() << " " << fPart[0].Y() << " " << fPart[0].Z() << " " << fPart[0].T() << endl;
            cout << "TLVec 1: " << fPart[1].X() << " " << fPart[1].Y() << " " << fPart[1].Z() << " " << fPart[1].T() << endl;
            cout << "TLVec 2: " << fPart[2].X() << " " << fPart[2].Y() << " " << fPart[2].Z() << " " << fPart[2].T() << endl;
	  }

   }//event loop

  // fFile->Close();

  TCanvas* can2 = new TCanvas("can2","Jake Data Plots",0,0,2000,1000);
  can2->Divide(2,2);

  can2->cd(1);
  wP[0]->Draw();
  can2->cd(2);
  wP[1]->Draw();
  can2->cd(3);
  wP[2]->Draw();
  can2->cd(4);
  w23->Draw();

  can2->Print("JakePlots.pdf"); 
  can2->Print("JakePlots.eps");
}
