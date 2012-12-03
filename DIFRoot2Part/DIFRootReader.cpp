#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <utility>
#include "DIFRootReader.hpp"
#include "TParticle.h"

DIFRootReader::DIFRootReader(const std::string inConfigFile, const bool binned=false):fBinned(binned){

  fFile = new TFile(inConfigFile.c_str());
  fTree = (TTree*) fFile->Get("data");
  fParticles = new TClonesArray("TParticle");
  fTree->GetBranch("Particles")->SetAutoDelete(false);
  fTree->SetBranchAddress("Particles",&fParticles);
  fFile->cd();

  fmaxEvents=fTree->GetEntries();
  fEvent=0;

  if(fBinned)
    bin();
}

DIFRootReader::~DIFRootReader(){
  fFile->Close();
  delete fParticles;
  delete fFile;
  //delete _myFcn;
}

const int DIFRootReader::getEvent(const int i, PWAEvent& inEvent){

  if(i>=0) {fEvent=i;}
  else {fEvent++;}

  fParticles->Clear();
  fTree->GetEntry(fEvent);

  // Get number of particle in TClonesrray
  unsigned int nParts = fParticles->GetEntriesFast();
  if(nParts!=2) return 0;

  TParticle* part1 = (TParticle*) fParticles->At(0); //pip
  TParticle* part2 = (TParticle*) fParticles->At(1); //pim
  if(!part1 || !part2) return 0;
  TLorentzVector in1, in2;
  part1->Momentum(in1);
  part2->Momentum(in2);
  std::vector<PWAParticle> out;
  out.push_back(PWAParticle(in1.X(), in1.Y(), in1.Z(), in1.E()));
  out.push_back(PWAParticle(in2.X(), in2.Y(), in2.Z(), in2.E()));

  //shared_ptr<PWAEvent> tmp(new PWAEvent());
  //inEvent = make_shared<PWAEvent>();
  for(unsigned int part=0; part<out.size(); part++)
    inEvent.addParticle(out.at(part));

  return 1;
}

const int DIFRootReader::getBin(const int i, double& m12, double& weight){
  if(!fBinned) return -1;

  m12 = fBins[i].first;
  weight = fBins[i].second;

  return 1;
}

const int DIFRootReader::getEvent(const int i,TLorentzVector& in1, TLorentzVector& in2, double& inm12){

  //TLorentzVector pPip,pPim,pPm12;
  //TLorentzVector V(0,0,0,0);
  //double m12sq;

  //fFile->cd();
  //TRandom3 rando;
  if(i>=0) {fEvent=i;}
  else {fEvent++;}
  fTree->GetEntry(fEvent);
  // Get number of particle in TClonesrray
  unsigned int nParts = fParticles->GetEntriesFast();
  if(nParts!=2) return -1;

  TParticle* part1 = (TParticle*) fParticles->At(0); //pip
  TParticle* part2 = (TParticle*) fParticles->At(1); //pim
  if(!part1 || !part2) return 0;
  part1->Momentum(in1);
  part2->Momentum(in2);

  //cout << part2->Energy() << endl;
  inm12=(in1+in2).Mag2();

  return 1;
}

void DIFRootReader::bin(){
  double min=-1, max=-1;
  fmaxBins=500; //TODO setter, consturctor
  TLorentzVector in1, in2;

  //initialize min & max
  fTree->GetEntry(0);
  unsigned int nParts = fParticles->GetEntriesFast();
  if(nParts!=2) return;
  TParticle* part1 = (TParticle*) fParticles->At(0); //pip
  TParticle* part2 = (TParticle*) fParticles->At(1); //pim
  if(!part1 || !part2) return;
  part1->Momentum(in1);
  part2->Momentum(in2);
  double inm12=(in1+in2).Mag2();
  min = max = inm12;

  //find min and max
  for(unsigned int evt=1; evt<fmaxEvents; evt++){
    fTree->GetEntry(evt);
    unsigned int nParts = fParticles->GetEntriesFast();
    if(nParts!=2) return;

    part1 = (TParticle*) fParticles->At(0); //pip
    part2 = (TParticle*) fParticles->At(1); //pim
    if(!part1 || !part2) return;
    part1->Momentum(in1);
    part2->Momentum(in2);

    inm12=(in1+in2).Mag2();

    if( min>inm12 ) min=inm12;
    if( max<inm12 ) max=inm12;
  }

  //initialize bins with weight zero
  min = sqrt(min);
  max = sqrt(max);
  double step = (max-min)/(double)fmaxBins;
  for(unsigned int bin=0; bin<fmaxBins; bin++){
    fBins[bin] = std::make_pair(min+bin*step,0.0);
  }

  //fill bins
  for(unsigned int evt=0; evt<fmaxBins; evt++){
    fTree->GetEntry(evt);
    unsigned int nParts = fParticles->GetEntriesFast();
    if(nParts!=2) return;

    TParticle* part1 = (TParticle*) fParticles->At(0); //pip
    TParticle* part2 = (TParticle*) fParticles->At(1); //pim
    if(!part1 || !part2) return;
    part1->Momentum(in1);
    part2->Momentum(in2);
    inm12=sqrt((in1+in2).Mag2());

    int bin = (int)((inm12-min)/step);
    fBins[bin].second += 1.;
  }

}
