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

//Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

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
  std::string file="test/3Part-4vecs.root";
  AmplitudeSetup ini("test/JPSI_ypipi.xml");
  std::cout << "Load Modules" << std::endl;
  std::shared_ptr<Data> myReader(new RootReader(file, false,"data"));
  std::shared_ptr<Data> myPHSPReader(new RootReader(file, false,"mc"));
  std::shared_ptr<Amplitude> amps(new AmpSumIntensity(M, Br, m1, m2, m3, ini));
  std::shared_ptr<Estimator> esti(new MinLogLH(amps, myReader, myPHSPReader));
  std::shared_ptr<Optimizer> opti(new MinuitIF(esti));

  // Initiate parameters
  ParameterList par;
  amps->fillStartParVec(par); //perfect startvalues
  std::cout << "LH mit optimalen intensitäten: " << esti->controlParameter(par) << std::endl;
  for(unsigned int i=0; i<par.GetNDouble(); i++){
    par.GetDoubleParameter(i).SetValue(300./(double)(i+1));
    par.GetDoubleParameter(i).SetError(10.);
  }
  std::cout << "LH mit folgenden intensitäten: " << esti->controlParameter(par) << std::endl;

  for(unsigned int i=0; i<par.GetNDouble(); i++){
      std::cout << "Parameter " << i << " = " << par.GetDoubleParameter(i).GetValue() << std::endl;
  }

 // std::cout << "Fixing 5 of 7 parameters " << std::endl;
  //for(unsigned int i=2; i<par.GetNDouble(); i++){
  //    par.GetDoubleParameter(i).FixParameter(true);
  //  }

  std::cout << "Start Fit" << std::endl;
  double genResult = opti->exec(par);
  std::cout << "Final LH = " << genResult << std::endl;

  //Plot result
  TH2D* bw12 = new TH2D("bw12","inv. mass-sq of particles 1&2 FitResult",1000,0.,10.,1000,0.,10.);
  bw12->GetXaxis()->SetTitle("m_{12}^{2} / GeV^{2}");
  bw12->GetXaxis()->CenterTitle();
  bw12->GetYaxis()->SetTitle("#");
  bw12->GetYaxis()->CenterTitle();
  TH2D* bw13 = new TH2D("bw13","inv. mass-sq of particles 1&3 FitResult",1000,0.,10.,1000,0.,10.);
  bw13->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
  bw13->GetXaxis()->CenterTitle();
  bw13->GetYaxis()->SetTitle("#");
  bw13->GetYaxis()->CenterTitle();
  TH2D* bw23 = new TH2D("bw23","inv. mass-sq of particles 2&3 FitResult",1000,0.,10.,1000,0.,10.);
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
      bw12FIT->Fill(masssq12,masssq13,amps->intensity(x,par));
      bw13FIT->Fill(masssq13,masssq12,amps->intensity(x,par));
      bw23FIT->Fill(masssq23,masssq12,amps->intensity(x,par));
  }

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
  output.Write();
  output.Close();

  return 0;
}
