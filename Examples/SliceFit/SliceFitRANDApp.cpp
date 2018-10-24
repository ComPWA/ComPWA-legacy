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
//! Test-Application for full fit with simple BW-dalitz-model.
/*!
 * @file DalitzFitApp.cpp
 * This tiny application tests a dalitz-fit procedure with a simple resonance
 * model. It uses the simple LH-estimator MinLogLH, it reads data
 * via the root-reader module RootReader and uses the intensity provided by
 * the Breit-Wigner-Sum  physics module AmplitudeSum. The optimization of the
 * parameters is done with the Minuit2 module MinuitIF. As result the
 * optimized parameters are printed to the terminal.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

// Root header files go here
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

#include "Data/RootReader/RootReader.hpp"
//Core header files go here
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"

// ComPWA header files go here
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Estimator/SliceFitUB/SliceFitUB.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"

// Stefans header files go here
#include "ROOTPlotHelper.hpp"

const Double_t M = 3.096916; // GeV/c² (J/psi+)
const Double_t Br = 0.000093; // GeV/c² (width)
const Double_t m1 = 0.; // GeV/c² (gamma)
const Double_t m2 = 0.139570; // GeV/c² (pi)
const Double_t m3 = 0.139570; // GeV/c² (pi)
//const Double_t c = 299792458.; // m/s
const Double_t PI = 3.14159; // m/s

unsigned int nFitEvents=100000-1;
unsigned int nStartEvent=0;
unsigned int nBins=160;
unsigned int nF0=3;
unsigned int nF2=2;

unsigned int testSlice = 15, nRand = 20;

using namespace ComPWA;
using namespace ComPWA::Physics::AmplitudeSum;
using Physics::DPKinematics::DalitzKinematics;

using Physics::AmplitudeSum::AmpSumIntensity;
using DataReader::RootReader::RootReader;
//using Estimator::SliceFitUB::SliceFitUB;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  boost::log::core::get()->set_filter(trivial::severity >= trivial::debug); //setting log level
  BOOST_LOG_TRIVIAL(info) << "  ComPWA Copyright (C) 2013  Mathias Michel ";
  BOOST_LOG_TRIVIAL(info) << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt";
  BOOST_LOG_TRIVIAL(info) << std::endl;

  DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(DalitzKinematics::createInstance("J/psi","gamma","pi0","pi0"));
	//DPKinematics kin("J/psi","gamma","pi0","pi0");
	//DPKinematics kin("D0","gamma","K-","K+");
	//static dataPoint* point = dataPoint::instance(kin);

  bool resultGen = true;

  std::string file="test/3Part-4vecs_1M_PDG.root";

	const char* pPath = getenv("COMPWA_DIR");
	std::string path = "";
	try{
	  path = std::string(pPath);
	}catch(std::logic_error& ex){
	  BOOST_LOG_TRIVIAL(error)<<"Environment Variable COMPWA_DIR not set?"<<std::endl;
	}
	std::string resoFile=path+"/test/JPSI_ypipi.xml";
    boost::property_tree::ptree pt;
    read_xml(resoFile, pt, boost::property_tree::xml_parser::trim_whitespace);

  BOOST_LOG_TRIVIAL(info)<< "Load Modules";
  std::shared_ptr<RootReader> myReader(new RootReader(file, "data"));
  std::shared_ptr<RootReader> myPHSPReader(new RootReader(file, "mc"));
  std::shared_ptr<AmpSumIntensity> amps(
		  new AmpSumIntensity(
				  "amp",
				  normStyle::none,
				  std::shared_ptr<Efficiency>(new UnitEfficiency()), nFitEvents)
  );
  amps->Configure(pt);

  // Initiate parameters
  ParameterList par;
  std::shared_ptr<Estimator::SliceFitUB::SliceFitUB> esti;
  amps->FillParameterList(par); //perfect startvalues
  esti = std::static_pointer_cast<Estimator::SliceFitUB::SliceFitUB>(
      Estimator::SliceFitUB::SliceFitUB::createInstance(
				  amps, myReader, myPHSPReader, par, nStartEvent,
				  nFitEvents, nBins, nF0, nF2
		  )
  );


  //unsigned int nSlices = nBins-(nBins/20.);
  ParameterList slicePars;
  //for(unsigned int i=0; i<nSlices; i++){
 //   std::string sliceName = "S"+std::to_string(i);
    std::shared_ptr<DoubleParameter> tmpA = std::shared_ptr<DoubleParameter>(new DoubleParameter("P0",0.75,0.1,100.));
    std::shared_ptr<DoubleParameter> tmpB = std::shared_ptr<DoubleParameter>(new DoubleParameter("P1",1.,-15,15));
    std::shared_ptr<DoubleParameter> tmpC = std::shared_ptr<DoubleParameter>(new DoubleParameter("P2",0.,0,15));
    std::shared_ptr<DoubleParameter> tmpD = std::shared_ptr<DoubleParameter>(new DoubleParameter("P3",0.1,-15,15));
    std::shared_ptr<DoubleParameter> tmpE = std::shared_ptr<DoubleParameter>(new DoubleParameter("P4",1.,0,50));
    //std::shared_ptr<DoubleParameter> tmpF = std::shared_ptr<DoubleParameter>(new DoubleParameter("P5",1.,-10,10));
    tmpA->FixParameter(true);
    //tmpC->FixParameter(true);
    //tmpE->FixParameter(true);
    //tmpB->FixParameter(true);
    //tmpD->FixParameter(true);
    tmpB->SetError(1.);
    tmpC->SetError(1.);
    tmpD->SetError(1.);
    tmpE->SetError(1.);
    //tmpA->UseBounds(true);
    //tmpB->UseBounds(true);
    //tmpC->UseBounds(true);
    //tmpD->UseBounds(true);
    //tmpE->UseBounds(true);
    slicePars.AddParameter(tmpA);
    slicePars.AddParameter(tmpB);
    slicePars.AddParameter(tmpC);
    slicePars.AddParameter(tmpD);
    slicePars.AddParameter(tmpE);
    //slicePars.AddParameter(tmpF);


    TH2D* phspA = new TH2D("phspTOT","phspTOT",100,0,10,100,0,10);
    TH2D* phspB = new TH2D("phspSpin0","phspSpin0",100,0,10,100,0,10);
    TH2D* phspD = new TH2D("phspSpin2","phspSpin2",100,0,10,100,0,10);
    TH2D* phspC = new TH2D("phspAll","phspAll",100,0,10,100,0,10);
    std::complex<double> reso[2];
    reso[0]=std::complex<double>(10,0);
    reso[1]=std::complex<double>(0,0);
    std::complex<double> resoC[2];
    resoC[0]=std::complex<double>(0,0);
    resoC[1]=std::complex<double>(20,0);
    std::complex<double> resoTOT[2];
    resoTOT[0]=std::complex<double>(0,0);
    resoTOT[1]=std::complex<double>(0,0);
    std::complex<double> resoFull[2];
    resoFull[0]=std::complex<double>(22,0);
    resoFull[1]=std::complex<double>(15,0);
    //std::cout << " " << reso[0] << "   " << reso[1] << std::endl;
    for(unsigned int i=0; i<95; i++){
      for(unsigned int j=0; j<95; j++){
      dataPoint point;
      point.setVal("m23sq",i/10.); point.setVal("m13sq",j/10.);
      //std::cout << " " << amps->sliceIntensity(point, par, reso, 2) << "   " << amps->sliceIntensity(point, par, resoTOT, 2) << std::endl;
      phspA->SetBinContent(i,j,amps->sliceIntensity(point, par, resoTOT, 2,1.,nF0,nF2));
      phspB->SetBinContent(i,j,amps->sliceIntensity(point, par, reso, 2,1.,nF0,nF2));
      phspC->SetBinContent(i,j,amps->sliceIntensity(point, par, resoFull, 2,1.,nF0,nF2));
      phspD->SetBinContent(i,j,amps->sliceIntensity(point, par, resoC, 2,1.,nF0,nF2));
    }}


    //slicePars.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("P6",1.)));
  //}

    std::shared_ptr<Optimizer::Optimizer> opti(new Optimizer::Minuit2::MinuitIF(esti, slicePars));

    esti->setSlice(33);


  BOOST_LOG_TRIVIAL(info) << "LH with start parameters slice 33: " << esti->controlParameter(slicePars);

  /*double startInt[par.GetNDouble()], optiInt[par.GetNDouble()];
  for(unsigned int i=0; i<par.GetNDouble(); i++){
    std::shared_ptr<DoubleParameter> tmp = par.GetDoubleParameter(i);
    optiInt[i] = tmp->GetValue();
    if(i<0 || i>9 || i%2==1){ //omega's and f0 fixed
   // if(i<2 || i>3 || i%2==1){ //omega's and f0 fixed
      tmp->FixParameter(true);
    }else{
      tmp->SetValue(tmp->GetValue()/((i+1)));
      tmp->SetError(std::shared_ptr<ParError<double>>(new SymError<double>(tmp->GetValue())));
      if(!tmp->GetValue()) tmp->SetError(std::shared_ptr<ParError<double>>(new SymError<double>(1.)));
    }
    startInt[i] = tmp->GetValue();
  }*/

 // std::cout << "Fixing 5 of 7 parameters " << std::endl;
  //for(unsigned int i=2; i<par.GetNDouble(); i++){
  //    par.GetDoubleParameter(i).FixParameter(true);
  //  }

  BOOST_LOG_TRIVIAL(info) << "Start multiple Fit of slice " << testSlice << " with " << nRand <<" random start parameter";
  std::vector<std::complex<double> > p1,p2, e1, e2;
  std::vector<double> invMass, norm, norme, lh;
  std::vector<std::shared_ptr<TH1D> > histData, histModel, histModelCl;
  //std::vector<std::shared_ptr<TH1D> > histData, histModel, histModelCl;
  TRandom3 r;
  for(int i=0; i<nRand; i++){


    double tmpMass = esti->setSlice(testSlice);
    BOOST_LOG_TRIVIAL(debug) << "InvMass Slice " << testSlice << " " << tmpMass ;
    std::shared_ptr<FitResult> genResult = opti->exec(slicePars);
    genResult->writeText("FitResultsAllSlicesUB_RAND.txt");

    histData.push_back(esti->getSliceHist());
    histModel.push_back(esti->getAmpSlHist());
    histModelCl.push_back(esti->getAmpClHist());

    histData.at(histData.size()-1)->GetXaxis()->SetTitle("m^{2}(#pi^{0} #gamma) [GeV/c^{2}]^{2}");
    histModel.at(histModel.size()-1)->GetXaxis()->SetTitle("m^{2}(#pi^{0} #gamma) [GeV/c^{2}]^{2}");
    histModelCl.at(histModelCl.size()-1)->GetXaxis()->SetTitle("m^{2}(#pi^{0} #gamma) [GeV/c^{2}]^{2}");
    histData.at(histData.size()-1)->GetYaxis()->SetTitle("Entries");
    histModel.at(histModel.size()-1)->GetYaxis()->SetTitle("Entries");
    histModelCl.at(histModelCl.size()-1)->GetYaxis()->SetTitle("Entries");
    histData.at(histData.size()-1)->SetStats(0);
    histModel.at(histModel.size()-1)->SetStats(0);
    histModelCl.at(histModelCl.size()-1)->SetStats(0);

    double parN(slicePars.GetDoubleParameter(0)->GetValue());
    double parE(slicePars.GetDoubleParameter(0)->GetError());
    std::complex<double> parCA(slicePars.GetDoubleParameter(1)->GetValue(), slicePars.GetDoubleParameter(2)->GetValue());
    std::complex<double> parCB(slicePars.GetDoubleParameter(3)->GetValue(), slicePars.GetDoubleParameter(4)->GetValue());
    std::complex<double> parEA(slicePars.GetDoubleParameter(1)->GetError(), slicePars.GetDoubleParameter(2)->GetError());
    std::complex<double> parEB(slicePars.GetDoubleParameter(3)->GetError(), slicePars.GetDoubleParameter(4)->GetError());
    p1.push_back(parCA);
    p2.push_back(parCB);
    e1.push_back(parEA);
    e2.push_back(parEB);
    invMass.push_back(i);
    norm.push_back(std::fabs(parN));
    norme.push_back(parE);
    lh.push_back(genResult->getResult());
    //p3.push_back(slicePars.GetDoubleParameter(3)->GetValue());
    //p4.push_back(slicePars.GetDoubleParameter(4)->GetValue());

    //if(i==30) i = 10000;
    genResult.reset();

    BOOST_LOG_TRIVIAL(info) << "RandFit " << i << " reset par" ;
    for(unsigned int i=0; i<slicePars.GetNDouble(); i++){
      std::shared_ptr<DoubleParameter> tmp = slicePars.GetDoubleParameter(i);
      if(!tmp->IsFixed()){
        tmp->SetValue(r.Uniform(0,10));
        tmp->SetError(tmp->GetValue());
      }
    }
  }

  BOOST_LOG_TRIVIAL(debug) << "Results";
  for(unsigned int i=0; i<p1.size(); i++){
    BOOST_LOG_TRIVIAL(debug) << "Rand " << i << " " << std::abs(p1[i]) << " " << std::abs(p2[i]);
  }

  if(!resultGen) return 0;

  //Plot result
  TGraphErrors par_N(nRand);
  TGraphErrors par_r0(nRand);
  TGraphErrors par_r2(nRand);
  TGraphErrors par_p0(nRand);
  TGraphErrors par_p2(nRand);

  TGraphErrors par_x0(nRand);
  TGraphErrors par_y0(nRand);
  TGraphErrors par_x2(nRand);
  TGraphErrors par_y2(nRand);

  TGraphErrors par_xy0(nRand);
  TGraphErrors par_xy2(nRand);

  TGraphErrors par_LH(nRand);

 // NeatPlotting::PlotBundle bundlexy0;

  //std::vector<NeatPlotting::PlotLabel> labelsxy0;
  for(unsigned int i=0; i<p1.size(); i++){
    double phi0=std::arg(p1[i]), phi2=std::arg(p2[i]);
    double phi0e, phi2e, abs0e, abs2e;
    while(phi0<0){
      phi0+=2*3.14159;
    };
    while(phi2<0){
      phi2+=2*3.14159;
    };

    double xx, yy, xexe, yeye;
    xx = p1[i].real()*p1[i].real(); yy = p1[i].imag()*p1[i].imag(); xexe = e1[i].real()*e1[i].real(); yeye = e1[i].imag()*e1[i].imag();
    phi0e = sqrt((xx*yeye + yy*xexe)/(xx+yy));
    abs0e = sqrt((xx*xexe + yy*yeye)/(xx+yy));
    xx = p2[i].real()*p2[i].real(); yy = p2[i].imag()*p2[i].imag(); xexe = e2[i].real()*e2[i].real(); yeye = e2[i].imag()*e2[i].imag();
    phi2e = sqrt((xx*yeye + yy*xexe)/(xx+yy));
    abs2e = sqrt((xx*xexe + yy*yeye)/(xx+yy));

    BOOST_LOG_TRIVIAL(debug) << "Slice " << i+(nBins/40.) << " " << std::abs(p1[i]) << " " << std::abs(p2[i]);
    BOOST_LOG_TRIVIAL(debug) << "Slice " << i+(nBins/40.) << " " << abs0e << " " << abs2e;

    par_r0.SetPoint(i,invMass[i],std::abs(p1[i]));
    par_p0.SetPoint(i,invMass[i],phi0);
    par_r0.SetPointError(i,0,abs0e);
    par_p0.SetPointError(i,0,phi0e);

    par_x0.SetPoint(i,invMass[i],p1[i].real());
    par_y0.SetPoint(i,invMass[i],p1[i].imag());
    par_x0.SetPointError(i,0,e1[i].real());
    par_y0.SetPointError(i,0,e1[i].imag());

    par_r2.SetPoint(i,invMass[i],std::abs(p2[i]));
    par_p2.SetPoint(i,invMass[i],phi2);
    par_r2.SetPointError(i,0,abs2e);
    par_p2.SetPointError(i,0,phi2e);

    par_x2.SetPoint(i,invMass[i],p2[i].real());
    par_y2.SetPoint(i,invMass[i],p2[i].imag());
    par_x2.SetPointError(i,0,e2[i].real());
    par_y2.SetPointError(i,0,e2[i].imag());

    par_xy0.SetPoint(i,p1[i].real(),p1[i].imag());
    par_xy2.SetPoint(i,p2[i].real(),p2[i].imag());
    par_xy0.SetPointError(i,e1[i].real(),e1[i].imag());
    par_xy2.SetPointError(i,e2[i].real(),e2[i].imag());
    //NeatPlotting::PlotLabel tmplabelxy0(std::to_string(i),NeatPlotting::TextStyle());
    //tmplabelxy0.setAbsolutionPosition(p1[i].real()+0.2,p1[i].imag()+0.2);
   // labelsxy0.push_back(tmplabelxy0);

    par_LH.SetPoint(i,invMass[i],lh[i]);

    par_N.SetPoint(i,invMass[i],norm[i]);
    par_N.SetPointError(i,0,norme[i]);
  }

  //bundlexy0.addGraph(par_xy0, NeatPlotting::DataObjectStyle());
 // for(unsigned int i=0; i<labelsxy0.size(); i++){
  //  bundlexy0.plot_decoration.labels.push_back(labelsxy0[i]);
  //}

  //NeatPlotting::PlotStyle myStyle;
  //bundlexy0.plot_axis.x_axis_title = "Real Spin 0";
 // bundlexy0.plot_axis.y_axis_title = "Imag Spin 0";

  //TCanvas bla("Argand Spin 0");
  //bundlexy0.drawOnCurrentPad(myStyle);


  par_r0.Draw();
  par_r2.Draw();
  par_p0.Draw();
  par_p2.Draw();
  par_r0.SetTitle("");
  par_r2.SetTitle("");
  par_p0.SetTitle("");
  par_p2.SetTitle("");
  par_r0.GetYaxis()->SetTitle("Magnitude Spin 0");
  par_p0.GetYaxis()->SetTitle("Phase Spin 0 /rad");
  par_r0.GetXaxis()->SetTitle("Run ID");
  par_p0.GetXaxis()->SetTitle("Run ID");
  par_r2.GetYaxis()->SetTitle("Magnitude Spin 2");
  par_p2.GetYaxis()->SetTitle("Phase Spin 2 /rad");
  par_r2.GetXaxis()->SetTitle("Run ID");
  par_p2.GetXaxis()->SetTitle("Run ID");

  par_x0.Draw();
  par_x2.Draw();
  par_y0.Draw();
  par_y2.Draw();
  par_x0.SetTitle("");
  par_x2.SetTitle("");
  par_y0.SetTitle("");
  par_y2.SetTitle("");
  par_x0.GetYaxis()->SetTitle("Real Spin 0");
  par_y0.GetYaxis()->SetTitle("Imag Spin 0");
  par_x0.GetXaxis()->SetTitle("Run ID");
  par_y0.GetXaxis()->SetTitle("Run ID");
  par_x2.GetYaxis()->SetTitle("Real Spin 2");
  par_y2.GetYaxis()->SetTitle("Imag Spin 2 /rad");
  par_x2.GetXaxis()->SetTitle("Run ID");
  par_y2.GetXaxis()->SetTitle("Run ID");

  par_xy0.Draw();
  par_xy2.Draw();
  par_xy0.SetTitle("");
  par_xy2.SetTitle("");
  par_xy0.GetXaxis()->SetTitle("Real Spin 0");
  par_xy0.GetYaxis()->SetTitle("Imag Spin 0");
  par_xy2.GetXaxis()->SetTitle("Real Spin 2");
  par_xy2.GetYaxis()->SetTitle("Imag Spin 2");

  par_LH.Draw();
  par_LH.SetTitle("");
  par_LH.GetXaxis()->SetTitle("RUN ID");
  par_LH.GetYaxis()->SetTitle("LH");

  TFile output("test/FitResultJPSISLICEUB_RAND.root","RECREATE","ROOT_Tree");
  par_N.Write("Constant");
  par_r0.Write("Magn Spin 0");
  par_r2.Write("Magn Spin 2");
  par_p0.Write("Arg Spin 0");
  par_p2.Write("Arg Spin 2");
  par_x0.Write("Real Spin 0");
  par_x2.Write("Real Spin 2");
  par_y0.Write("Imag Spin 0");
  par_y2.Write("Imag Spin 2");
  //phspA->Write();
 // phspB->Write();
 // phspC->Write();
 // phspD->Write();
  par_xy0.Write("Argand Spin 0");
  par_xy2.Write("Argand Spin 2");
  par_LH.Write("LH");
  output.mkdir("RandFit");
  output.cd("RandFit");
  for(unsigned int h=0; h<histData.size(); h++){
    histData[h]->Write();
    histModel[h]->Write();
    histModelCl[h]->Write();
  }
  output.Write();
  output.Close();

  BOOST_LOG_TRIVIAL(info) << "Done";

  return 0;
}
