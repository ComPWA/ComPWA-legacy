// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Root-Headers
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TParticlePDG.h"

#include "Core/Kinematics.hpp"
#include "Core/Generator.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "DataReader/RootReader/RootReader.hpp"

namespace ComPWA {
namespace DataReader {

RootReader::RootReader() {}

RootReader::RootReader(TTree *tr, int size) : Data(false) { read(tr, size); }

RootReader::RootReader(const std::string inRootFile,
                       const std::string inTreeName, int size)
    : Data(false) {

  TFile *fFile = new TFile(inRootFile.c_str());

  if (fFile->IsZombie())
    throw std::runtime_error("RootReader::RootReader() | "
                             "Can't open data file: " +
                             inRootFile);

  TTree *fTree = (TTree *)fFile->Get(inTreeName.c_str());

  if (!fTree)
    throw std::runtime_error("RootReader::RootReader() | Tree \"" + inTreeName +
                             "\" can not be opened from file " + inRootFile +
                             "! ");
  read(fTree, size);

  fFile->Close();
}

RootReader::~RootReader() {}

RootReader *RootReader::clone() const { return new RootReader(*this); }

RootReader *RootReader::emptyClone() const { return new RootReader(); }

void RootReader::read(TTree *fTree, double readSize) {

  // TTree branch variables
  TClonesArray *fParticles;
  double feventWeight;
  double feventEff;
  int fCharge;
  int fFlavour;

  fParticles = new TClonesArray("TParticle");
  fTree->GetBranch("Particles")->SetAutoDelete(false);
  fTree->SetBranchAddress("Particles", &fParticles);
  fTree->SetBranchAddress("eff", &feventEff);
  fTree->SetBranchAddress("weight", &feventWeight);
  fTree->SetBranchAddress("charge", &fCharge);
  fTree->SetBranchAddress("flavour", &fFlavour);

  if (readSize <= 0 || readSize > fTree->GetEntries())
    readSize = fTree->GetEntries();

  for (unsigned int evt = 0; evt < readSize; evt++) {
    Event tmp;
    fParticles->Clear();
    fTree->GetEntry(evt);

    // Get number of particle in TClonesrray
    unsigned int nParts = fParticles->GetEntriesFast();

    TParticle *partN;
    TLorentzVector inN;
    for (unsigned int part = 0; part < nParts; part++) {
      partN = 0;
      partN = (TParticle *)fParticles->At(part);
      if (!partN)
        continue;
      partN->Momentum(inN);
      int charge = partN->GetPDG()->Charge();
      if (charge != 0)
        charge /= std::fabs(charge);
      tmp.addParticle(Particle(inN.X(), inN.Y(), inN.Z(), inN.E(),
                               partN->GetPdgCode(), charge));
    } // particle loop
    tmp.setWeight(feventWeight);
    tmp.setCharge(fCharge);
    tmp.setEfficiency(feventEff);

    Events.push_back(tmp);
    if (feventWeight > MaximumWeight)
      MaximumWeight = feventWeight;
  } // end event loop
}

void RootReader::writeData(std::string fileName, std::string treeName) {

  LOG(INFO) << "RootReader::writeData() | Writing current "
               "vector of events to file "
            << fileName;

  TFile *fFile = new TFile(fileName.c_str(), "RECREATE");
  if (fFile->IsZombie())
    throw std::runtime_error(
        "RootReader::writeData() | Can't open data file: " + fileName);

  // TTree branch variables
  TClonesArray *fParticles;
  double feventWeight;
  double feventEff;
  int fCharge;
  int fFlavour;

  TTree *fTree = new TTree(treeName.c_str(), treeName.c_str());
  unsigned int numPart = Events.at(0).numParticles();
  fParticles = new TClonesArray("TParticle", numPart);
  fTree->Branch("Particles", &fParticles);
  fTree->Branch("weight", &feventWeight, "weight/D");
  fTree->Branch("eff", &feventEff, "weight/D");
  fTree->Branch("charge", &fCharge, "charge/I");
  fTree->Branch("flavour", &fFlavour, "flavour/I");
  TClonesArray &partArray = *fParticles;

  auto it = Events.begin();
  for (; it != Events.end(); ++it) {
    fParticles->Clear();
    feventWeight = (*it).weight();
    fCharge = (*it).charge();
    feventEff = (*it).efficiency();

    TLorentzVector motherMomentum(0, 0, 0, (*it).cmsEnergy());
    for (unsigned int i = 0; i < numPart; i++) {
      Particle oldParticle = (*it).particle(i);
      TLorentzVector oldMomentum(oldParticle.px(), oldParticle.py(),
                                 oldParticle.pz(), oldParticle.e());
      new (partArray[i]) TParticle(oldParticle.pid(), 1, 0, 0, 0, 0,
                                   oldMomentum, motherMomentum);
    }
    fTree->Fill();
  }
  fTree->Write("", TObject::kOverwrite, 0);
  fFile->Close();
  return;
}

} // namespace DataReader
} // namespace ComPWA
