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
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
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

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt" << std::endl;
  std::cout << std::endl;

  DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(DalitzKinematics::createInstance("J/psi","gamma","pi0","pi0"));
	//DPKinematics kin("J/psi","gamma","pi0","pi0");
	//DPKinematics kin("D0","gamma","K-","K+");
	//static dataPoint* point = dataPoint::instance(kin);

  bool useFctTree = false;

  std::string file="test/3Part-4vecs.root";

	const char* pPath = getenv("COMPWA_DIR");
	std::string path = "";
	try{
	  path = std::string(pPath);
	}catch(std::logic_error){
	  BOOST_LOG_TRIVIAL(error)<<"Environment Variable COMPWA_DIR not set?"<<std::endl;
	}
	std::string resoFile=path+"/test/JPSI_ypipi.xml";
	AmplitudeSetup ini(resoFile);

  std::cout << "Load Modules" << std::endl;
  std::shared_ptr<Data> myReader(new RootReader(file, false,"data"));
  std::shared_ptr<Data> myPHSPReader(new RootReader(file, false,"mc"));
  std::shared_ptr<Amplitude> amps(new AmpSumIntensity(ini,std::shared_ptr<Efficiency>(new UnitEfficiency()), myReader->getNEvents()));

  //std::shared_ptr<Amplitude> amps(new AmpSumIntensity(M, Br, m1, m2, m3,"J/psi","gamma","pi0","pi0", ini));
  // Initiate parameters
  ParameterList par;
  std::shared_ptr<ControlParameter> esti;
  if(!useFctTree){//using tree?
    amps->fillStartParVec(par); //perfect startvalues
    //std::cout << "Pars: " << par.GetNDouble() << std::endl;
    esti = MinLogLH::createInstance(amps, myReader, myPHSPReader);
  }else{
    std::shared_ptr<FunctionTree> physicsTree = amps->functionTree(par);
    esti = MinLogLH::createInstance(physicsTree, myReader, myPHSPReader);
  }
  //std::shared_ptr<ControlParameter> esti = MinLogLH::createInstance(amps, myReader, myPHSPReader);
  //std::shared_ptr<Estimator> esti(new MinLogLH(amps, myReader, myPHSPReader));
  std::shared_ptr<Optimizer> opti(new MinuitIF(esti, par));

  ParameterList test;
  if(useFctTree)
    if(!amps->functionTree(test))
      return 1;

  //return 0;

  std::cout << "LH with optimal parameters: " << esti->controlParameter(par) << std::endl;
  double startInt[par.GetNDouble()], optiInt[par.GetNDouble()];
  for(unsigned int i=0; i<par.GetNDouble(); i++){
    std::shared_ptr<DoubleParameter> tmp = par.GetDoubleParameter(i);
    optiInt[i] = tmp->GetValue();
    if(i<2 || i>3){ //omega's and f0 fixed
      tmp->FixParameter(true);
    }else{
      tmp->SetValue(tmp->GetValue()/((i+1)));
      tmp->SetError(std::shared_ptr<ParError<double>>(new SymError<double>(tmp->GetValue())));
      if(!tmp->GetValue()) tmp->SetError(std::shared_ptr<ParError<double>>(new SymError<double>(1.)));
    }
    startInt[i] = tmp->GetValue();
  }
  std::cout << "LH with following parameters: " << esti->controlParameter(par) << std::endl;



  for(unsigned int i=0; i<par.GetNDouble(); i++){
      std::cout << par.GetDoubleParameter(i)->GetName() << " = " << par.GetDoubleParameter(i)->GetValue() << std::endl;
  }

 // std::cout << "Fixing 5 of 7 parameters " << std::endl;
  //for(unsigned int i=2; i<par.GetNDouble(); i++){
  //    par.GetDoubleParameter(i).FixParameter(true);
  //  }

  std::cout << "Start Fit" << std::endl;
  std::shared_ptr<FitResult> genResult = opti->exec(par);
  std::cout << "Final LH = " << genResult->getResult() << std::endl;

  std::cout << "Optimierte intensitäten: " << esti->controlParameter(par) << std::endl;

  for(unsigned int i=0; i<par.GetNDouble(); i++){
      std::cout << par.GetDoubleParameter(i)->GetName() << " = " << par.GetDoubleParameter(i)->GetValue();
      std::cout << "   [ start: " << startInt[i] << " ,";
      std::cout << " optimal: " << optiInt[i] << " ]" << std::endl;
  }
/*
  //Plot result
  TH2D* bw12 = new TH2D("bw12","inv. mass-sq of particles 1&2 Generated",1000,0.,10.,1000,0.,10.);
  bw12->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
  bw12->GetXaxis()->CenterTitle();
  bw12->GetYaxis()->SetTitle("#");
  bw12->GetYaxis()->CenterTitle();
  TH2D* bw13 = new TH2D("bw13","inv. mass-sq of particles 1&3 Generated",1000,0.,10.,1000,0.,10.);
  bw13->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
  bw13->GetXaxis()->CenterTitle();
  bw13->GetYaxis()->SetTitle("#");
  bw13->GetYaxis()->CenterTitle();
  TH2D* bw23 = new TH2D("bw23","inv. mass-sq of particles 2&3 Generated",1000,0.,10.,1000,0.,10.);
  bw23->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
  bw23->GetXaxis()->CenterTitle();
  bw23->GetYaxis()->SetTitle("#");
  bw23->GetYaxis()->CenterTitle();

  RootReader myReaderPHSP(file, false,"mc");
  unsigned int maxEventsPHSP = myReaderPHSP.getNEvents();
  //double masssq12PHSP, masssq13PHSP, masssq23PHSP;
  TH2D* bw12PHSP = new TH2D("bw12PHSP","inv. mass-sq of particles 1&2 PHSP",1000,0.,10.,1000,0.,10.);
  bw12PHSP->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
  bw12PHSP->GetXaxis()->CenterTitle();
  bw12PHSP->GetYaxis()->SetTitle("#");
  bw12PHSP->GetYaxis()->CenterTitle();
  TH2D* bw13PHSP = new TH2D("bw13PHSP","inv. mass-sq of particles 1&3 PHSP",1000,0.,10.,1000,0.,10.);
  bw13PHSP->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
  bw13PHSP->GetXaxis()->CenterTitle();
  bw13PHSP->GetYaxis()->SetTitle("#");
  bw13PHSP->GetYaxis()->CenterTitle();
  TH2D* bw23PHSP = new TH2D("bw23PHSP","inv. mass-sq of particles 2&3 PHSP",1000,0.,10.,1000,0.,10.);
  bw23PHSP->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
  bw23PHSP->GetXaxis()->CenterTitle();
  bw23PHSP->GetYaxis()->SetTitle("#");
  bw23PHSP->GetYaxis()->CenterTitle();

  TH2D* bw12FIT = new TH2D("bw12FIT","inv. mass-sq of particles 1&2 FitResult",1000,0.,10.,1000,0.,10.);
  bw12FIT->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
  bw12FIT->GetXaxis()->CenterTitle();
  bw12FIT->GetYaxis()->SetTitle("#");
  bw12FIT->GetYaxis()->CenterTitle();
  TH2D* bw13FIT = new TH2D("bw13FIT","inv. mass-sq of particles 1&3 FitResult",1000,0.,10.,1000,0.,10.);
  bw13FIT->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
  bw13FIT->GetXaxis()->CenterTitle();
  bw13FIT->GetYaxis()->SetTitle("#");
  bw13PHSP->GetYaxis()->CenterTitle();
  TH2D* bw23FIT = new TH2D("bw23FIT","inv. mass-sq of particles 2&3 FitResult",1000,0.,10.,1000,0.,10.);
  bw23FIT->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
  bw23FIT->GetXaxis()->CenterTitle();
  bw23FIT->GetYaxis()->SetTitle("#");
  bw23FIT->GetYaxis()->CenterTitle();

  TH2D* bw12DIFF = new TH2D("bw12DIFF","inv. mass-sq of particles 1&2 FitResult",1000,0.,10.,1000,0.,10.);
  bw12DIFF->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
  bw12DIFF->GetXaxis()->CenterTitle();
  bw12DIFF->GetYaxis()->SetTitle("#");
  bw12DIFF->GetYaxis()->CenterTitle();
  TH2D* bw13DIFF = new TH2D("bw13DIFF","inv. mass-sq of particles 1&3 FitResult",1000,0.,10.,1000,0.,10.);
  bw13DIFF->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
  bw13DIFF->GetXaxis()->CenterTitle();
  bw13DIFF->GetYaxis()->SetTitle("#");
  bw13PHSP->GetYaxis()->CenterTitle();
  TH2D* bw23DIFF = new TH2D("bw23DIFF","inv. mass-sq of particles 2&3 FitResult",1000,0.,10.,1000,0.,10.);
  bw23DIFF->GetXaxis()->SetTitle("m_{23}^{2} / GeV^{2}");
  bw23DIFF->GetXaxis()->CenterTitle();
  bw23DIFF->GetYaxis()->SetTitle("#");
  bw23DIFF->GetYaxis()->CenterTitle();

  double masssq12, masssq13, masssq23;
  for(unsigned int i = 0; i < myReader->getNEvents(); i++){
      Event event(myReader->getEvent(i));

      //myReader.getEvent(-1, a, b, masssq);
      //if(!myReader.getEvent(i, event)) continue; TODO: try exception
      if(!event.getNParticles() == 3) continue;
      //if(!event) continue;
      //cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles() << endl;
      const Particle &a(event.getParticle(0));
      const Particle &b(event.getParticle(1));
      const Particle &c(event.getParticle(2));
      masssq12 = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);
      masssq13 = pow(a.E+c.E,2) - pow(a.px+c.px ,2) - pow(a.py+c.py ,2) - pow(a.pz+c.pz ,2);
      masssq23 = pow(b.E+c.E,2) - pow(b.px+c.px ,2) - pow(b.py+c.py ,2) - pow(b.pz+c.pz ,2);

      bw12->Fill(masssq12,masssq13);
      bw13->Fill(masssq13,masssq12);
      bw23->Fill(masssq23,masssq12);
  }

 // TH2D* bw12DIFF = (TH2D*)bw12->Clone("bw12DIFF");
 // TH2D* bw13DIFF = (TH2D*)bw13->Clone("bw13DIFF");
 // TH2D* bw23DIFF = (TH2D*)bw23->Clone("bw23DIFF");

  for(unsigned int i = 0; i < maxEventsPHSP; i++){
      Event event(myReaderPHSP.getEvent(i));

      //myReader.getEvent(-1, a, b, masssq);
      //if(!myReader.getEvent(i, event)) continue; TODO: try exception
      if(!event.getNParticles() == 3) continue;
      //if(!event) continue;
      //cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles() << endl;
      const Particle &a(event.getParticle(0));
      const Particle &b(event.getParticle(1));
      const Particle &c(event.getParticle(2));
      masssq12 = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);
      masssq13 = pow(a.E+c.E,2) - pow(a.px+c.px ,2) - pow(a.py+c.py ,2) - pow(a.pz+c.pz ,2);
      masssq23 = pow(b.E+c.E,2) - pow(b.px+c.px ,2) - pow(b.py+c.py ,2) - pow(b.pz+c.pz ,2);

      bw12PHSP->Fill(masssq12,masssq13);
      bw13PHSP->Fill(masssq13,masssq12);
      bw23PHSP->Fill(masssq23,masssq12);

      vector<double> x;
      x.push_back(sqrt(masssq23));
      x.push_back(sqrt(masssq13));
      x.push_back(sqrt(masssq12));
      //bw12FIT->Fill(masssq12,masssq13,1000*amps->intensity(x,par));
    //  bw13FIT->Fill(masssq13,masssq12,1000*amps->intensity(x,par));
     // bw23FIT->Fill(masssq23,masssq12,1000*amps->intensity(x,par));
  }


  //Generation
  TRandom3 rando;
  TLorentzVector W(0.0, 0.0, 0.0, M);//= beam + target;

  //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses[3] = { m1, m2, m2} ;
  TGenPhaseSpace event;
  event.SetDecay(W, 3, masses);

  TLorentzVector *pGamma,*pPip,*pPim,pPm23,pPm13,pPm12;
  double weight, m23sq, m13sq, m12sq, maxTest=0;
  ParameterList paras(par);
  if(useFctTree){
    paras.RemoveDouble("ma"); paras.RemoveDouble("mb"); paras.RemoveDouble("mc");
  }

  cout << "Einschwingen" << endl;
  for(unsigned int schwing=0; schwing<10*myReader->getNEvents(); schwing++){
      weight = event.Generate();

      pGamma = event.GetDecay(0);
      pPip    = event.GetDecay(1);
      pPim    = event.GetDecay(2);

	  pPm23 = *pPim + *pPip;
	  pPm13 = *pGamma + *pPim;
	  pPm12 = *pGamma + *pPip;

      m23sq=pPm23.M2(); m13sq=pPm13.M2(); m12sq=pPm12.M2();

      //		m12sq = kin.getThirdVariableSq(m23sq,m13sq);
      		point->setMsq(3,m12sq); point->setMsq(4,m13sq); point->setMsq(5,m23sq);
      //		m12sq=M*M+m1*m1+m2*m2+m3*m3-m13sq-m23sq;
      		if( abs(m12sq-kin.getThirdVariableSq(m23sq,m13sq))>0.01 ){
      			std::cout<<m12sq<<" "<<kin.getThirdVariableSq(m23sq,m13sq)<<std::endl;
      			std::cout<<"   " <<m23sq<<" "<<m13sq<<" "<<m12sq<<std::endl;
      		}

      //call physics module
      vector<double> x;
      x.push_back(sqrt(m23sq));
      x.push_back(sqrt(m13sq));
      x.push_back(sqrt(m12sq));
      ParameterList intensL = amps->intensity(x, paras);
      double AMPpdf = intensL.GetDoubleParameter(0)->GetValue();
      //double AMPpdf = testBW.intensity(x, minPar);


      //mb.setVal(m13);
      //double m13pdf = totAmp13.getVal();//fun_combi2->Eval(m13);
      if(maxTest<(weight*AMPpdf))
        maxTest=(weight*AMPpdf);

  }

  maxTest*=1.1;
  cout << "Start generation of y pi0 pi0 Dalitz Result" << endl;
  unsigned int i = 0;
  do{
        weight = event.Generate();

        pGamma = event.GetDecay(0);
        pPip    = event.GetDecay(1);
        pPim    = event.GetDecay(2);

        pPm23 = *pPim + *pPip;
        pPm13 = *pGamma + *pPim;
        pPm12 = *pGamma + *pPip;

        m23sq=pPm23.M2(); m13sq=pPm13.M2(); m12sq=pPm12.M2();

        point->setMsq(3,m12sq); point->setMsq(4,m13sq); point->setMsq(5,m23sq);

       // m12sq=M*M-m13sq-m23sq;
        //if(m12sq<0){
          //cout << tmpm12_sq << "\t" << M*M << "\t" << m13_sq << "\t" << m23_sq << endl;
          //continue;
        //  m12sq=0.0001;
       // }

        TParticle fparticleGam(22,1,0,0,0,0,*pGamma,W);
        TParticle fparticlePip(211,1,0,0,0,0,*pPip,W);
        TParticle fparticlePim(-211,1,0,0,0,0,*pPim,W);

        //call physics module
	dataPoint dataP; dataP.setVal("m23sq",m23sq);	dataP.setVal("m13sq",m13sq);
//        vector<double> x;
//        x.push_back(sqrt(m23sq));
//        x.push_back(sqrt(m13sq));
//        x.push_back(sqrt(m12sq));
        ParameterList intensL = amps->intensity(dataP, paras);
        double AMPpdf = intensL.GetDoubleParameter(0)->GetValue();
        //double AMPpdf = amps->intensity(x, par);

        double test = rando.Uniform(0,maxTest);

        //mb.setVal(m13);
        //double m13pdf = totAmp13.getVal();//fun_combi2->Eval(m13);
        if(maxTest<(weight*AMPpdf))
          cout << "Einschwingen zu kurz!" << endl;
        if(test<(weight*AMPpdf)){
          i++;

          bw12FIT->Fill(m12sq,m13sq);
          bw13FIT->Fill(m13sq,m12sq);
          bw23FIT->Fill(m23sq,m12sq);
        }

    }while(i<100000);

  bw12DIFF->Add(bw12,1.0);
  bw23DIFF->Add(bw23,1.0);
  bw13DIFF->Add(bw13,1.0);
  bw12DIFF->Divide(bw12FIT);
  bw23DIFF->Divide(bw23FIT);
  bw13DIFF->Divide(bw13FIT);
 // bw12DIFF = new TH2D(*bw12 - *bw12FIT);
 // bw23DIFF = new TH2D(*bw23 - *bw23FIT);
 // bw13DIFF = new TH2D(*bw13 - *bw13FIT);

  TFile output("test/FitResultJPSI.root","RECREATE","ROOT_Tree");
  bw12->Write();
  bw13->Write();
  bw23->Write();
  bw12PHSP->Write();
  bw13PHSP->Write();
  bw23PHSP->Write();
  bw12FIT->Write();
  bw13FIT->Write();
  bw23FIT->Write();
  bw12DIFF->Write();
  bw13DIFF->Write();
  bw23DIFF->Write();
  output.Write();
  output.Close();
*/
  cout << "Done" << endl;

  return 0;
}
