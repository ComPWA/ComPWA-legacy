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
#include "Core/Kinematics.hpp"
#include "Core/Generator.hpp"
#include "TParticle.h"
#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
using namespace boost::log;

void RootReader::Clear(){
	fEvents.clear();
}
std::shared_ptr<Data> RootReader::rndSubSet(unsigned int size, std::shared_ptr<Generator> gen){
	std::shared_ptr<Data> newSample(new RootReader(fileName, true,"test",false));
	unsigned int totalSize = getNEvents();
	unsigned int newSize = totalSize;
	/* 1th method: new Sample has exact size, but possibly events are added twice.
	 * We would have to store all used events in a vector and search the vector at every event -> slow
	 */
	//unsigned int t=0;
	//unsigned int d=0;
	//while(t<newSize){
	//	d = (unsigned int) gen->getUniform()*totalSize;
	//	newSample->pushEvent(fEvents[d]);
	//	t++;
	//}

	/* 2nd method: events are added once only, but total size of sample varies with sqrt(N) */
	for(unsigned int i=0; i<totalSize; i++){//count how many events are not within PHSP
		dataPoint point(fEvents[i]);
		if(!Kinematics::instance()->isWithinPhsp(point)) newSize--;
	}
	double threshold = (double)size/newSize; //calculate threshold
	for(unsigned int i=0; i<totalSize; i++){
		dataPoint point(fEvents[i]);
		if(!Kinematics::instance()->isWithinPhsp(point)) continue;
		if(gen->getUniform()<threshold)
			newSample->pushEvent(fEvents[i]);
	}
	return newSample;
}

void RootReader::read(){
	fParticles = new TClonesArray("TParticle");
	fTree->GetBranch("Particles")->SetAutoDelete(false);
	fTree->SetBranchAddress("Particles",&fParticles);
	fTree->SetBranchAddress("eff",&feventEff);
	fTree->SetBranchAddress("weight",&feventWeight);
	fTree->SetBranchAddress("charge",&fCharge);
	fTree->SetBranchAddress("flavour",&fFlavour);
	//	fEvent=0;
	bin();
	storeEvents();

}
RootReader::RootReader(TTree* tr, const bool binned=false) : fBinned(binned){
	fTree = tr;
	fFile = 0; //need to do this to avoid seg. violation when destructor is called
	read();
}
RootReader::RootReader(const std::string inRootFile, const bool binned,
		const std::string inTreeName, const bool readFlag)
:fBinned(binned),_readFlag(readFlag),fileName(inRootFile),treeName(inTreeName){
	//	fEvent=0;
	if(!readFlag) return;
	fFile = new TFile(fileName.c_str());
	fTree = (TTree*) fFile->Get(treeName.c_str());
	read();
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

//void RootReader::writeToFile(){
//	if(_readFlag){
//		std::cout<<"RootReader: trying to write, but RootReader is marked as readonly! Dont write!"<<std::endl;
//		return;
//	}
//	fFile = new TFile(fileName.c_str(),"RECREATE");
//	fTree = new TTree(treeName.c_str(),treeName.c_str());
//	TParticle* part = 0;
//	fTree->Branch("Particles","Particles",&part,64000,0);
//	//loop
//	for(int i=0; i<=fEvents.size();i++){
//
//
//	}
//	return;
//}
const std::vector<std::string>& RootReader::getVariableNames(){
	if(!fVarNames.size()){ //TODO: init
		fVarNames.push_back("dataname1");
		fVarNames.push_back("dataname2");
	}
	return fVarNames;
}

Event& RootReader::getEvent(const int i){
	//Event outEvent;

	//	if(i>=0) {fEvent=i;}
	//	else {fEvent++;}

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

allMasses RootReader::getMasses(){
	if(!fEvents.size()) return allMasses();
	unsigned int nParts = fEvents.at(0).getNParticles();
	//	BOOST_LOG_TRIVIAL(debug)<<"RootReader::getMasses() #particles: "<<nParts;

	//determine invMass combinations
	unsigned int nMasses=0;
	std::vector<std::pair<unsigned int, unsigned int> > ids;
	for(unsigned int i=0; i<nParts; i++)
		for(unsigned int j=i+1; j<nParts; j++){
			nMasses++;
			ids.push_back(std::make_pair(i+1,j+1));
		}
	//	BOOST_LOG_TRIVIAL(debug)<<"RootReader::getMasses() #invMasses: "<<nMasses;

	unsigned int nSkipped =0; //count events which are outside PHSP boundary
	unsigned int nFilled=0; //count events which are outside PHSP boundary

	//	allMasses result(nMasses, fEvents.size(), ids);
	allMasses result(nMasses, ids);
	//calc and store inv masses
	for(unsigned int evt=0; evt<fEvents.size(); evt++){
		if( result.Fill(fEvents.at(evt)) ) nFilled++;
		else nSkipped++;

		//		Event tmp = fEvents.at(evt);
		//		dataPoint point(tmp);
		//
		//		// Check number of particle in TClonesrray and if event is within PHSP boundary
		//		if( nParts != tmp.getNParticles()  || !Kinematics::instance()->isWithinPhsp(point)){
		//			nSkipped++;
		////			result.nEvents--;
		//			continue;
		//		}
		////		result.eff.at(evt) = tmp.getEfficiency();
		////		result.weight.at(evt) = tmp.getWeight();
		//		result.eff.push_back(tmp.getEfficiency());
		//		result.weight.push_back(tmp.getWeight());
		//
		//		for(unsigned int pa=0; pa<nParts; pa++){
		//			for(unsigned int pb=pa+1; pb<nParts; pb++){
		//				const Particle &inA(tmp.getParticle(pa));
		//				const Particle &inB(tmp.getParticle(pb));
		//				double mymass_sq = inA.invariantMass(inB);
		////				(result.masses_sq.at(std::make_pair(pa+1,pb+1))).at(evt) = mymass_sq;
		//				(result.masses_sq.at(std::make_pair(pa+1,pb+1))).push_back(mymass_sq);
		//			}//particle loop B
		//		}//particle loop A
		//		result.nEvents++;
		//		nFilled++;

	}//event loop
	if(nSkipped)
		BOOST_LOG_TRIVIAL(debug)<<"RootReader::getMasses() "<<nSkipped
		<<"("<<(double)nSkipped/fEvents.size()*100<<"%) data points are "
				"outside the PHSP boundary. We skip these points!";

	//	std::cout<<"after      "<<result.masses_sq.at(std::make_pair(2,3)).size()<<std::endl;
	return result;
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

	for(unsigned int evt=0; evt<fTree->GetEntries(); evt++){
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
			tmp.addParticle(Particle(inN.X(), inN.Y(), inN.Z(), inN.E(),partN->GetPdgCode()));
		}//particle loop
		tmp.setWeight(feventWeight);
		tmp.setCharge(fCharge);
		tmp.setFlavour(fFlavour);
		tmp.setEfficiency(feventEff);

		fEvents.push_back(tmp);
	}//event loop

}

void RootReader::writeData(){
	BOOST_LOG_TRIVIAL(info) << "RootReader: writing current vector of events to file "<<fileName;
	TFile* ff = new TFile(fileName.c_str(),"UPDATE");
	fTree = new TTree(treeName.c_str(),treeName.c_str());
	unsigned int numPart = fEvents[0].getNParticles();
	fParticles = new TClonesArray("TParticle",numPart);
	fTree->Branch("Particles",&fParticles);
	fTree->Branch("weight",&feventWeight,"weight/D");
	fTree->Branch("eff",&feventEff,"weight/D");
	fTree->Branch("charge",&fCharge,"charge/I");
	fTree->Branch("flavour",&fFlavour,"flavour/I");
	TClonesArray &partArray = *fParticles;

	TLorentzVector motherMomentum(0,0,0,Kinematics::instance()->getMotherMass());
	for(std::vector<Event>::iterator it=fEvents.begin(); it!=fEvents.end(); it++){
		fParticles->Clear();
		feventWeight = (*it).getWeight();
		fCharge= (*it).getCharge();
		fFlavour= (*it).getFlavour();
		feventEff = (*it).getEfficiency();
		for(unsigned int i=0; i<numPart; i++){
			const Particle oldParticle = (*it).getParticle(i);
			TLorentzVector oldMomentum(oldParticle.px,oldParticle.py,oldParticle.pz,oldParticle.E);
			new(partArray[i]) TParticle(oldParticle.pid,1,0,0,0,0,oldMomentum,motherMomentum);
		}
		fTree->Fill();
	}
	fTree->Write("",TObject::kOverwrite,0);
	ff->Close();
	return;
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
	for(unsigned int evt=1; evt<fTree->GetEntries(); evt++){
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
