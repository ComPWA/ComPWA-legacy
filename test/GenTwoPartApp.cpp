//! Application to generate a two-particle intensity.
/*!
 * @file GenTwoPartApp.cpp
 * This application uses the simple Breit-Wigner Physics-Module and the root
 * phase-space generator to generate a file with two-particle events.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

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
#include "Physics/BreitWigner/BreitWigner.hpp"
#include "PWAParameter.hpp"
#include "PWAGenericPar.hpp"

using namespace std;

const unsigned int MaxEvents = 100000;

//constants
const Double_t M = 3.096916; // GeV/c² (J/Psi mass)
//const Double_t B = 0.1851; // GeV/c² (f2 width)
const Double_t m1 = 0.139570; // GeV/c² (pi)
const Double_t PI = 3.14159;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  unsigned int i=0;
  TRandom3 rando;

  //Simple Breit-Wigner Physics-Module setup
  shared_ptr<BreitWigner> testBW(new BreitWigner(0.,5.));
  vector<shared_ptr<PWAParameter> > minPar;
  minPar.push_back(shared_ptr<PWAGenericPar<double> >(new PWAGenericPar<double>(1.5,0.5,2.5,0.1)));
  minPar.push_back(shared_ptr<PWAGenericPar<double> >(new PWAGenericPar<double>(0.3,0.1,0.2,0.01)));

  //Output File setup
  TFile output("test/2Part-4vecs.root","recreate");
  output.SetCompressionLevel(1); //try level 2 also
  TTree fTree ("data","Dalitz-Gen");
  TClonesArray *fEvt = new TClonesArray("TParticle");
  TClonesArray &ar = *fEvt;
  fTree.Branch("Particles",&fEvt);

  //Generation
  TLorentzVector *pPip,*pPim, pPm;
  double weight;
  cout << "Start generation of 2Particle Breit-Wigner" << endl;
  do{
      TLorentzVector W(0.0, 0.0, 0.0, rando.Uniform(2*m1,2.4));//= beam + target;

     //(Momentum, Energy units are Gev/C, GeV)
      Double_t masses[2] = { m1, m1} ;

      TGenPhaseSpace event;
      event.SetDecay(W, 2, masses);

      weight = event.Generate();

      pPip    = event.GetDecay(0);
      pPim    = event.GetDecay(1);

      pPm = *pPim + *pPip;

      TParticle fparticlePip(211,1,0,0,0,0,*pPip,W);
      TParticle fparticlePim(-211,1,0,0,0,0,*pPim,W);

      //call physics module
      vector<double> x;
      x.push_back(pPm.Mag());
      double BWpdf = testBW->intensity(x, minPar); //TMath::BreitWigner(pPm.Mag(),1.2,0.2);

      double test = rando.Uniform(0,10);

      if(test<(weight*BWpdf)){
        ar.Clear();
        new(ar[0]) TParticle(fparticlePip);
        new(ar[1]) TParticle(fparticlePim);
        //if((i/1000.-TMath::Floor(i/1000.))==0.) std::cout << "Energy: " << W.E() << std::endl;
        fTree.Fill();
        i++;
      }
  }while(i<MaxEvents);

  fTree.Print();
  fTree.Write();
  output.Close();

  cout << "Done ..." << endl << endl;

  return 0;
}
