//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <utility>
#include "DataReader/RootReader/RootReader.hpp"
#include "TParticle.h"

void RootReader::init(){
	fParticles = new TClonesArray("TParticle");
	fTree->GetBranch("Particles")->SetAutoDelete(false);
	fTree->SetBranchAddress("Particles",&fParticles);
	fmaxEvents=fTree->GetEntries();
	fEvent=0;
	bin();
	storeEvents();

}
RootReader::RootReader(TTree* tr, const bool binned=false) : fBinned(binned){
	fTree = tr;
	fFile = 0; //need to do this to avoid seg. violation when destructor is called
	init();
}
RootReader::RootReader(const std::string inRootFile, const bool binned,
		const std::string inTreeName, const bool readFlag)
:fBinned(binned),_readFlag(readFlag),fileName(inRootFile),treeName(inTreeName){
std::cout<<"sdfdsf"<<std::endl;
	fEvent=0;
	if(!readFlag) return;
	fFile = new TFile(fileName.c_str());
	fTree = (TTree*) fFile->Get(treeName.c_str());
	init();
	//	fParticles = new TClonesArray("TParticle");
	//	fTree->GetBranch("Particles")->SetAutoDelete(false);
	//	fTree->SetBranchAddress("Particles",&fParticles);
	//	fFile->cd();
	//
	//	fmaxEvents=fTree->GetEntries();
	//	fEvent=0;
	//
	//	//if(fBinned)
	//	bin();
	//	storeEvents();

	fFile->Close();
}
RootReader::~RootReader(){
	//fFile->Close();
//	delete fParticles;
//	delete fFile;
	//delete _myFcn;
}

void RootReader::writeToFile(){
	if(_readFlag){
		std::cout<<"RootReader: trying to write, but RootReader is marked as readonly! Exit!"<<std::endl;
		exit(1);
	}
	fFile = new TFile(fileName.c_str(),"RECREATE");
	fTree = new TTree(treeName.c_str(),treeName.c_str());
	TParticle* part = 0;
	fTree->Branch("Particles","Particles",&part,64000,0);

	//loop
	for(int i=0; i<=fEvents.size();i++){


	}

	return;
}
const std::vector<std::string>& RootReader::getVariableNames(){
	if(!fVarNames.size()){ //TODO: init
		fVarNames.push_back("dataname1");
		fVarNames.push_back("dataname2");
	}
	return fVarNames;
}

const Event& RootReader::getEvent(const int i){
	//Event outEvent;

	if(i>=0) {fEvent=i;}
	else {fEvent++;}

	return fEvents.at(i);

	/*fParticles->Clear();
  fTree->GetEntry(fEvent);

  // Get number of particle in TClonesrray
  unsigned int nParts = fParticles->GetEntriesFast();

  TParticle* partN;
  TLorentzVector inN;
  for(unsigned int part=0; part<nParts; part++){
    partN = 0;
    partN = (TParticle*) fParticles->At(part);
    if(!partN) continue;
    partN->Momentum(inN);
    outEvent.addParticle(Particle(inN.X(), inN.Y(), inN.Z(), inN.E()));
  }

  if(nParts!=2) return 0;

  TParticle* part1 = (TParticle*) fParticles->At(0); //pip
  TParticle* part2 = (TParticle*) fParticles->At(1); //pim
  if(!part1 || !part2) return 0;
  TLorentzVector in1, in2;
  part1->Momentum(in1);
  part2->Momentum(in2);
  std::vector<Particle> out;
  out.push_back(Particle(in1.X(), in1.Y(), in1.Z(), in1.E()));
  out.push_back(Particle(in2.X(), in2.Y(), in2.Z(), in2.E()));

  //shared_ptr<Event> tmp(new Event());
  //inEvent = make_shared<Event>();
  for(unsigned int part=0; part<out.size(); part++)
    inEvent.addParticle(out.at(part));*/

	//return outEvent;
}

const int RootReader::getBin(const int i, double& m12, double& weight){
	if(!fBinned) return -1;

	m12 = fBins[i].first;
	weight = fBins[i].second;

	return 1;
}

/*const int RootReader::getEvent(const int i,TLorentzVector& in1, TLorentzVector& in2, double& inm12){

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
}*/

void RootReader::storeEvents(){

	for(unsigned int evt=0; evt<fmaxEvents; evt++){
		Event tmp;
		fParticles->Clear();
		fTree->GetEntry(evt);

		// Get number of particle in TClonesrray
		unsigned int nParts = fParticles->GetEntriesFast();

		TParticle* partN;
		TLorentzVector inN;
		for(unsigned int part=0; part<nParts; part++){
			partN = 0;
			partN = (TParticle*) fParticles->At(part);
			if(!partN) continue;
			partN->Momentum(inN);
			tmp.addParticle(Particle(inN.X(), inN.Y(), inN.Z(), inN.E()));
		}//particle loop

		fEvents.push_back(tmp);
	}//event loop

}


void RootReader::bin(){
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
