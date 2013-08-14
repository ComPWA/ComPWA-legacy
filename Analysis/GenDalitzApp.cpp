//! Application to generate J/Psi -> g pi pi.
/*!
 * @file GenDalitzApp.cpp
 * This application uses the Breit-Wigner-Amplitude-Sum module and a
 * phase-space generator to generate a file with J/Psi -> g pi pi events.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

// Boost header files go here
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

// Root header files go here
#include "TMath.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TFile.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

// Physics Interface header files go here
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

#include "Analysis/PlotData.hpp"

using namespace std;

const unsigned int MaxEvents = 150000;

//constants
const Double_t M = 1.86486; // GeV/c² (D0)
const Double_t Br = 0.0; // GeV/c² (width)
//const Double_t m1 = 0.14; // GeV/c² (K_S^0)
//const Double_t m2 = 0.14; // GeV/c² (K^-)
//const Double_t m3 = 0.14; // GeV/c² (K^+)
const Double_t m1 = 0.497614; // GeV/c² (K_S^0)
const Double_t m2 = 0.493677; // GeV/c² (K^-)
const Double_t m3 = 0.493677; // GeV/c² (K^+)
const Double_t PI = 3.14159;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  unsigned int i=0, mc=0;
  TRandom3 rando;

  std::string outFile="gen-out.root";
  //load resonances
  AmplitudeSetup ini("Analysis/DKsKKRes.xml");
  cout << "loaded file " << ini.getFileName() << " with " << ini.getResonances().size() << " resonances!" << std::endl;
  //Simple Breit-Wigner Physics-Module setup
  AmpSumIntensity testBW(M, Br, m1, m2, m3, ini);
  testBW.printAmps();
  ParameterList minPar;
  testBW.fillStartParVec(minPar);

  //Output File setup
  TTree fTree ("data","Dalitz-Gen");
  TClonesArray *fEvt = new TClonesArray("TParticle");
  //TClonesArray &ar = *fEvt;
  fTree.Branch("Particles",&fEvt);

  TTree fTreePHSP ("mc","Dalitz-Gen-PHSP");
  TClonesArray *fEvtPHSP = new TClonesArray("TParticle");
  //TClonesArray &ar = *fEvt;
  fTreePHSP.Branch("Particles",&fEvtPHSP);

  //Generation
  TLorentzVector W(0.0, 0.0, 0.0, M);//= beam + target;

  //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses[3] = { m1, m2, m3} ;

  TGenPhaseSpace event;
  event.SetDecay(W, 3, masses);

  TLorentzVector *p0,*p1,*p2,pPm23,pPm13,pPm12;
  double weight, m23sq, m13sq, m12sq, maxTest=0;
   int scale = (int) MaxEvents/10;
  cout << "Start generation of K_S0 K- K+ Dalitzplot" << endl;
  do{
      weight = event.Generate();

      p0 = event.GetDecay(0);
      p1    = event.GetDecay(1);
      p2    = event.GetDecay(2);

      pPm23 = *p2 + *p1;
      pPm13 = *p0 + *p2;
      pPm12 = *p0 + *p1;

      m23sq=pPm23.M2(); m13sq=pPm13.M2(); m12sq=pPm12.M2();

      //m12sq=M*M-m13sq-m23sq;//WRONG! Sum(s_i)=M^2+Sum(m_i^2)
      //m12sq = M*M+m1*m1+m2*m2+m3*m3-m13sq-m23sq;

      if(m12sq<0){
        //cout << tmpm12_sq << "\t" << M*M << "\t" << m13_sq << "\t" << m23_sq << endl;
        //continue;
        m12sq=0.0001;
      }

      TParticle fparticleGam(22,1,0,0,0,0,*p0,W);
      TParticle fparticlePip(211,1,0,0,0,0,*p1,W);
      TParticle fparticlePim(-211,1,0,0,0,0,*p2,W);

      //call physics module
      vector<double> x;
      x.push_back(sqrt(m23sq));
      x.push_back(sqrt(m13sq));
//      x.push_back(sqrt(m12sq));//use only two parameters
      double AMPpdf = testBW.intensity(x, minPar);
//      double test = rando.Uniform(0,5);
      double test = rando.Uniform(0,testBW.getMaxVal()*2);

      if(maxTest<(weight*AMPpdf))
        maxTest=(weight*AMPpdf);
      if(i<MaxEvents && test<(weight*AMPpdf)){
    	  if( (i % scale) == 0) {cout<<"Progress [%]: "<<(i/scale)*10<<endl;}
        i++;
        new((*fEvt)[0]) TParticle(fparticleGam);
        new((*fEvt)[1]) TParticle(fparticlePip);
        new((*fEvt)[2]) TParticle(fparticlePim);

        fTree.Fill();
      }

      if(mc<MaxEvents && test<weight){
        mc++;
        new((*fEvtPHSP)[0]) TParticle(fparticleGam);
        new((*fEvtPHSP)[1]) TParticle(fparticlePip);
        new((*fEvtPHSP)[2]) TParticle(fparticlePim);

        fTreePHSP.Fill();
      }
  }while(i<MaxEvents || mc<MaxEvents);
  plotData plot(std::string("dalitz"),outFile, std::shared_ptr<Data>(new RootReader(&fTree,0)));  plot.plot();
  plotData plotPhsp(std::string("dalitz_phsp"),outFile, std::shared_ptr<Data>(new RootReader(&fTreePHSP,0)));  plotPhsp.plot();

  TFile output(outFile.c_str(),"update");
  output.SetCompressionLevel(1); //try level 2 also
  fTree.Print();
  fTree.Write("",TObject::kOverwrite,0);
  fTreePHSP.Write("",TObject::kOverwrite,0);
  output.Close();

  cout << "Done ... " << maxTest << endl << endl;

  return 0;
}
