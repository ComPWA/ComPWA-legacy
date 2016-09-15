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
#include <stdexcept>

#include "DataReader/RootReader/RootReader.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Generator.hpp"

#include "TParticle.h"

#include "TParticlePDG.h"
#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
using namespace boost::log;

namespace ComPWA {
namespace DataReader {
namespace RootReader {

RootReader::RootReader(){
	fFile = 0; //need to do this to avoid seg. violation when destructor is called
}

RootReader::RootReader(TTree* tr, int size, const bool binned) :
						readSize(size), Data(binned)
{
	fTree = tr;
	fFile = 0; //need to do this to avoid seg. violation when destructor is called
	read();
}

RootReader::RootReader(const std::string inRootFile,
		const std::string inTreeName, int size, const bool binned) :
							fileName(inRootFile),treeName(inTreeName),
							readSize(size),  Data(binned)
{
	fFile = new TFile(fileName.c_str());

	if(fFile->IsZombie())
		throw std::runtime_error("RootReader::RootReader() | "
				"Can't open data file: "+inRootFile);

	fTree = (TTree*) fFile->Get(treeName.c_str());

	if(!fTree)
		throw std::runtime_error("RootReader::RootReader() | Tree \""+treeName
				+"\" can not be opened from file "+inRootFile+"! ");

	read();

	fFile->Close();
}

RootReader::~RootReader()
{

}

RootReader* RootReader::Clone() const
{
	return new RootReader(*this);
}

RootReader* RootReader::EmptyClone() const
{
	return new RootReader();
}

void RootReader::read()
{
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


void RootReader::storeEvents()
{
	if(readSize<=0||readSize> fTree->GetEntries() )
		readSize = fTree->GetEntries();

	for(unsigned int evt=0; evt<readSize; evt++){
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
			int charge = partN->GetPDG()->Charge();
			if(charge != 0) charge /= std::fabs(charge);
			tmp.addParticle(
					Particle( inN.X(), inN.Y(), inN.Z(), inN.E(),
							partN->GetPdgCode(), charge )
			);
		}//particle loop
		tmp.setWeight(feventWeight);
		tmp.setCharge(fCharge);
		tmp.setFlavour(fFlavour);
		tmp.setEfficiency(feventEff);

		fEvents.push_back(tmp);
		if(feventWeight > maxWeight) maxWeight = feventWeight;
	}//event loop
}

void RootReader::writeData(std::string file, std::string trName)
{
	if(file!="") fileName=file;
	if(trName!="") treeName=trName;

	BOOST_LOG_TRIVIAL(info) << "RootReader::writeData() | Writing current "
			"vector of events to file "<<fileName;

	TFile* fFile = new TFile(fileName.c_str(),"UPDATE");
	if(fFile->IsZombie())
		throw std::runtime_error("RootReader::writeData() | "
				"Can't open data file: "+fileName);

	fTree = new TTree(treeName.c_str(),treeName.c_str());
	unsigned int numPart = fEvents[0].getNParticles();
	fParticles = new TClonesArray("TParticle",numPart);
	fTree->Branch("Particles",&fParticles);
	fTree->Branch("weight",&feventWeight,"weight/D");
	fTree->Branch("eff",&feventEff,"weight/D");
	fTree->Branch("charge",&fCharge,"charge/I");
	fTree->Branch("flavour",&fFlavour,"flavour/I");
	TClonesArray &partArray = *fParticles;

	TLorentzVector motherMomentum(0,0,0,Kinematics::instance()->GetMotherMass());
	auto it = fEvents.begin();
	for( ; it!=fEvents.end(); ++it){
		fParticles->Clear();
		feventWeight = (*it).getWeight();
		fCharge= (*it).getCharge();
		fFlavour= (*it).getFlavour();
		feventEff = (*it).getEfficiency();
		for(unsigned int i=0; i<numPart; i++){
			const Particle oldParticle = (*it).getParticle(i);
			TLorentzVector oldMomentum(
					oldParticle.px,oldParticle.py,oldParticle.pz,oldParticle.E
			);
			new(partArray[i]) TParticle(
					oldParticle.pid,1,0,0,0,0,oldMomentum,motherMomentum
			);
		}
		fTree->Fill();
	}
	fTree->Write("",TObject::kOverwrite,0);
	fFile->Close();
	return;
}

void RootReader::bin()
{
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

} /* namespace RootReader */
} /* namespace DataReader */
} /* namespace ComPWA */
