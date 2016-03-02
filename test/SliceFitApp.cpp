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
#include "TMath.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

//Core header files go here
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Estimator/SliceFit/SliceFit.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"

const Double_t M = 3.096916; // GeV/c² (J/psi+)
const Double_t Br = 0.000093; // GeV/c² (width)
const Double_t m1 = 0.; // GeV/c² (gamma)
const Double_t m2 = 0.139570; // GeV/c² (pi)
const Double_t m3 = 0.139570; // GeV/c² (pi)
//const Double_t c = 299792458.; // m/s
const Double_t PI = 3.14159; // m/s

unsigned int nFitEvents=100000-1;
unsigned int nStartEvent=0;
unsigned int nBins=200;

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

	std::string file="test/3Part-4vecs.root";

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
	auto fitAmpPtr = new AmpSumIntensity(normStyle::none,
			std::shared_ptr<Efficiency>(new UnitEfficiency()), nFitEvents);
	fitAmpPtr->Configure(pt);
	std::shared_ptr<AmpSumIntensity> amps( fitAmpPtr );

	BOOST_LOG_TRIVIAL(info)<< "Load Modules";
	std::shared_ptr<RootReader> myReader(new RootReader(file, "data"));
	std::shared_ptr<RootReader> myPHSPReader(new RootReader(file, "mc"));

	// Initiate parameters
	ParameterList par;
	std::shared_ptr<SliceFit> esti;
	amps->copyParameterList(par); //perfect startvalues
	esti = std::static_pointer_cast<SliceFit>(
			SliceFit::createInstance(
					amps, myReader, myPHSPReader, par, nStartEvent, nFitEvents
					)
	);


	//unsigned int nSlices = nBins-(nBins/20.);
	ParameterList slicePars;
	//for(unsigned int i=0; i<nSlices; i++){
	//   std::string sliceName = "S"+std::to_string(i);
	std::shared_ptr<DoubleParameter> tmpA = std::shared_ptr<DoubleParameter>(new DoubleParameter("P0",0.05,0,100));
	std::shared_ptr<DoubleParameter> tmpB = std::shared_ptr<DoubleParameter>(new DoubleParameter("P1",10.,0,20));
	std::shared_ptr<DoubleParameter> tmpC = std::shared_ptr<DoubleParameter>(new DoubleParameter("P2",0.,0,20));
	std::shared_ptr<DoubleParameter> tmpD = std::shared_ptr<DoubleParameter>(new DoubleParameter("P3",0.,0,30));
	std::shared_ptr<DoubleParameter> tmpE = std::shared_ptr<DoubleParameter>(new DoubleParameter("P4",0.,0,20));
	//std::shared_ptr<DoubleParameter> tmpF = std::shared_ptr<DoubleParameter>(new DoubleParameter("P5",1.,-10,10));
	tmpA->FixParameter(true);
	//tmpC->FixParameter(true);
	// tmpE->FixParameter(true);
	//tmpB->FixParameter(true);
	//tmpD->FixParameter(true);
	tmpB->SetError(1.);
	tmpC->SetError(1.);
	tmpD->SetError(1.);
	tmpE->SetError(1.);
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
	for(unsigned int i=0; i<100; i++){
		for(unsigned int j=0; j<100; j++){
			dataPoint point;
			point.setVal("m23sq",i/10.); point.setVal("m13sq",j/10.);
			//std::cout << " " << amps->sliceIntensity(point, par, reso, 2) << "   " << amps->sliceIntensity(point, par, resoTOT, 2) << std::endl;
			phspA->SetBinContent(i,j,amps->sliceIntensity(point, par, resoTOT, 2));
			phspB->SetBinContent(i,j,amps->sliceIntensity(point, par, reso, 2));
			phspC->SetBinContent(i,j,amps->sliceIntensity(point, par, resoFull, 2));
			phspD->SetBinContent(i,j,amps->sliceIntensity(point, par, resoC, 2));
		}}


	//slicePars.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("P6",1.)));
	//}

	std::shared_ptr<Optimizer> opti(new MinuitIF(esti, slicePars));

	esti->setSlice(33);

	BOOST_LOG_TRIVIAL(info) << "Chi2 with start parameters slice 33: " << esti->controlParameter(slicePars);

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

	BOOST_LOG_TRIVIAL(info) << "Start Fit of all slices";
	std::vector<std::complex<double> > p1,p2;
	std::vector<double> invMass, norm;
	std::vector<std::shared_ptr<TH1D> > histData, histModel, histModelCl;
	for(int i=0+(nBins/40.); i<nBins-(nBins/20.)-1; i++){
		//  {unsigned int i=50;

		BOOST_LOG_TRIVIAL(info) << "Slice " << i << " reset par" ;
		for(unsigned int i=0; i<slicePars.GetNDouble(); i++){
			std::shared_ptr<DoubleParameter> tmp = slicePars.GetDoubleParameter(i);
			if(!tmp->IsFixed()){
				tmp->SetValue(1.);
				tmp->SetError(tmp->GetValue());
			}
		}

		double tmpMass = esti->setSlice(i);
		BOOST_LOG_TRIVIAL(debug) << "InvMass Slice " << i << " " << tmpMass ;
		std::shared_ptr<FitResult> genResult = opti->exec(slicePars);

		histData.push_back(esti->getSliceHist());
		histModel.push_back(esti->getAmpSlHist());
		histModelCl.push_back(esti->getAmpClHist());

		double parN(slicePars.GetDoubleParameter(0)->GetValue());
		std::complex<double> parCA(slicePars.GetDoubleParameter(1)->GetValue(), slicePars.GetDoubleParameter(2)->GetValue());
		std::complex<double> parCB(slicePars.GetDoubleParameter(3)->GetValue(), slicePars.GetDoubleParameter(4)->GetValue());
		p1.push_back(parCA);
		p2.push_back(parCB);
		invMass.push_back(tmpMass);
		norm.push_back(std::fabs(parN));
		//p3.push_back(slicePars.GetDoubleParameter(3)->GetValue());
		//p4.push_back(slicePars.GetDoubleParameter(4)->GetValue());

	}

	BOOST_LOG_TRIVIAL(debug) << "Results";
	for(unsigned int i=0; i<p1.size(); i++){
		BOOST_LOG_TRIVIAL(debug) << "Slice " << i+(nBins/40.) << " " << std::abs(p1[i]) << " " << std::abs(p2[i]);
	}

	if(!resultGen) return 0;

	//Plot result
	TGraphErrors par_N(nBins-3*(nBins/40.)-1);
	TGraphErrors par_r0(nBins-3*(nBins/40.)-1);
	TGraphErrors par_r2(nBins-3*(nBins/40.)-1);
	TGraphErrors par_p0(nBins-3*(nBins/40.)-1);
	TGraphErrors par_p2(nBins-3*(nBins/40.)-1);


	for(unsigned int i=0; i<p1.size(); i++){
		double phi0=std::arg(p1[i]), phi2=std::arg(p2[i]);
		while(phi0<0){
			phi0+=2*3.14159;
		};
		while(phi2<0){
			phi2+=2*3.14159;
		};

		par_r0.SetPoint(i,invMass[i],std::abs(p1[i]));
		par_p0.SetPoint(i,invMass[i],phi0);

		par_r2.SetPoint(i,invMass[i],std::abs(p2[i]));
		par_p2.SetPoint(i,invMass[i],phi2);

		par_N.SetPoint(i,invMass[i],norm[i]);
	}

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
	par_r0.GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}/c^{2}");
	par_p0.GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}/c^{2}");
	par_r2.GetYaxis()->SetTitle("Magnitude Spin 2");
	par_p2.GetYaxis()->SetTitle("Phase Spin 2 /rad");
	par_r2.GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}/c^{2}");
	par_p2.GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}/c^{2}");


	TFile output("test/FitResultJPSISLICE.root","RECREATE","ROOT_Tree");
	par_N.Write();
	par_r0.Write();
	par_r2.Write();
	par_p0.Write();
	par_p2.Write();
	phspA->Write();
	phspB->Write();
	phspC->Write();
	phspD->Write();
	output.mkdir("Slices");
	output.cd("Slices");
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
