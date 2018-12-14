#include "TTree.h"
#include "TParticle.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TClonesArray.h"

///
/// Example ROOT macro to demonstrate how to generate a TTree that can be 
/// handled by RootReader.
///
void createRootTTree(){

  TFile* tf = new TFile("test.root","RECREATE");
  TTree* tr = new TTree("test","test");
  TClonesArray parts("TParticle",3);
  double weight, eff;
  int flavour, charge;
  tr->Branch("Particles",&parts);
  tr->Branch("weight",&weight);
  tr->Branch("eff",&eff);
  tr->Branch("charge",&charge);
  tr->Branch("flavour",&flavour);

  TLorentzVector p4_cms(0,0,0,3.77);
  double masses[3] = { 0.5, 0.5, 0.5 };
  TGenPhaseSpace* gen = new TGenPhaseSpace();
  gen->SetDecay(p4_cms, 3, masses);

  TRandom3 rand(0);
  for(int i=0; i<100; i++){
	eff = rand.Uniform(0,1);
	weight = rand.Uniform(0,1);
	parts = TClonesArray("TParticle",3);
	gen->Generate();

	TLorentzVector p4_1 = *gen->GetDecay(0);
	TLorentzVector p4_2 = *gen->GetDecay(1);
	TLorentzVector p4_3 = *gen->GetDecay(2);
	new(parts[0]) TParticle(310,1,0,0,0,0,p4_1,p4_cms);
	new(parts[1]) TParticle(-321,1,0,0,0,0,p4_2,p4_cms);
	new(parts[2]) TParticle(321,1,0,0,0,0,p4_3,p4_cms);
	tr->Fill();
  }
  tr->Write();
  tf->Close();
}
