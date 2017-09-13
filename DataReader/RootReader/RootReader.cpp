// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "TParticle.h"
#include "TLorentzVector.h"
#include "TParticlePDG.h"

#include "Core/Kinematics.hpp"
#include "Core/Generator.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "DataReader/RootReader/RootReader.hpp"

using namespace boost::log;

namespace ComPWA {
namespace DataReader {

RootReader::RootReader() {}

RootReader::RootReader(TTree *tr, int size) : Data(false) {
  read(tr, size);
}

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

RootReader *RootReader::Clone() const { return new RootReader(*this); }

RootReader *RootReader::EmptyClone() const { return new RootReader(); }

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
      tmp.AddParticle(Particle(inN.X(), inN.Y(), inN.Z(), inN.E(),
                               partN->GetPdgCode(), charge));
    } // particle loop
    tmp.SetWeight(feventWeight);
    tmp.SetCharge(fCharge);
    tmp.SetFlavour(fFlavour);
    tmp.SetEfficiency(feventEff);

    fEvents.push_back(tmp);
    if (feventWeight > maxWeight)
      maxWeight = feventWeight;
  } // end event loop
}

void RootReader::WriteData(std::string fileName, std::string treeName) {

  LOG(info) << "RootReader::writeData() | Writing current "
               "vector of events to file "
            << fileName;

  TFile *fFile = new TFile(fileName.c_str(), "RECREATE");
  if (fFile->IsZombie())
    throw std::runtime_error(
        "RootReader::WriteData() | Can't open data file: " + fileName);

  // TTree branch variables
  TClonesArray *fParticles;
  double feventWeight;
  double feventEff;
  int fCharge;
  int fFlavour;

  TTree* fTree = new TTree(treeName.c_str(), treeName.c_str());
  unsigned int numPart = fEvents.at(0).GetNParticles();
  fParticles = new TClonesArray("TParticle", numPart);
  fTree->Branch("Particles", &fParticles);
  fTree->Branch("weight", &feventWeight, "weight/D");
  fTree->Branch("eff", &feventEff, "weight/D");
  fTree->Branch("charge", &fCharge, "charge/I");
  fTree->Branch("flavour", &fFlavour, "flavour/I");
  TClonesArray &partArray = *fParticles;

  auto it = fEvents.begin();
  for (; it != fEvents.end(); ++it) {
    fParticles->Clear();
    feventWeight = (*it).GetWeight();
    fCharge = (*it).GetCharge();
    fFlavour = (*it).GetFlavour();
    feventEff = (*it).GetEfficiency();

    TLorentzVector motherMomentum(0, 0, 0, (*it).GetCMSEnergy());
    for (unsigned int i = 0; i < numPart; i++) {
      const Particle oldParticle = (*it).GetParticle(i);
      TLorentzVector oldMomentum(oldParticle.GetPx(), oldParticle.GetPy(),
                                 oldParticle.GetPz(), oldParticle.GetE());
      new (partArray[i]) TParticle(oldParticle.GetPid(), 1, 0, 0, 0, 0,
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
