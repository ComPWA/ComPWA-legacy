// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Root-Headers
#include "TClonesArray.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TTree.h"

#include "Core/Generator.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "RootReader.hpp"

namespace ComPWA {
namespace Data {

RootReader::RootReader(const std::string TreeName_, int NumberEventsToProcess_)
    : TreeName(TreeName_), NumberEventsToProcess(NumberEventsToProcess_) {}

RootReader::~RootReader() {}

// RootReader *RootReader::clone() const { return new RootReader(*this); }

// RootReader *RootReader::emptyClone() const { return new RootReader(); }

std::shared_ptr<ComPWA::Data::Data>
RootReader::readData(const std::string &InputFilePath) const {

  TFile File(InputFilePath.c_str());
  if (File.IsZombie())
    throw std::runtime_error("RootReader::RootReader() | "
                             "Can't open data file: " +
                             InputFilePath);

  TTree *fTree = (TTree *)File.Get(TreeName.c_str());

  if (!fTree)
    throw std::runtime_error("RootReader::RootReader() | Tree \"" + TreeName +
                             "\" can not be opened from file " + InputFilePath +
                             "! ");

  // TTree branch variables
  TClonesArray Particles("TParticle");
  TClonesArray *pParticles(&Particles);
  double feventWeight;
  double feventEff;
  int fCharge;
  int fFlavour;

  fTree->GetBranch("Particles")->SetAutoDelete(false);
  fTree->SetBranchAddress("Particles", &pParticles);
  fTree->SetBranchAddress("eff", &feventEff);
  fTree->SetBranchAddress("weight", &feventWeight);
  fTree->SetBranchAddress("charge", &fCharge);
  fTree->SetBranchAddress("flavour", &fFlavour);

  unsigned int NumberEventsToRead(NumberEventsToProcess);
  if (NumberEventsToProcess <= 0 || NumberEventsToProcess > fTree->GetEntries())
    NumberEventsToRead = fTree->GetEntries();

  std::shared_ptr<ComPWA::Data::Data> data(new ComPWA::Data::Data);
  std::vector<ComPWA::Event> Events;
  Events.reserve(NumberEventsToRead);

  for (unsigned int i = 0; i < NumberEventsToRead; ++i) {
    Event evt;
    Particles.Clear();
    fTree->GetEntry(i);

    // Get number of particle in TClonesrray
    unsigned int nParts = Particles.GetEntriesFast();

    TParticle *partN;
    TLorentzVector inN;
    for (unsigned int part = 0; part < nParts; part++) {
      partN = 0;
      partN = (TParticle *)Particles.At(part);
      if (!partN)
        continue;
      partN->Momentum(inN);
      int charge = partN->GetPDG()->Charge();
      if (charge != 0)
        charge /= std::fabs(charge);
      evt.addParticle(Particle(inN.X(), inN.Y(), inN.Z(), inN.E(),
                               partN->GetPdgCode(), charge));
    } // particle loop
    evt.setWeight(feventWeight);
    evt.setCharge(fCharge);
    evt.setEfficiency(feventEff);

    Events.push_back(evt);
  } // end event loop
  data->add(Events);
  File.Close();
  return data;
}

void RootReader::writeData(std::shared_ptr<ComPWA::Data::Data> Data,
                           const std::string &OutputFilePath) const {

  LOG(INFO) << "RootReader::writeData() | Writing current "
               "vector of events to file "
            << OutputFilePath;

  TFile File(OutputFilePath.c_str(), "RECREATE");
  if (File.IsZombie())
    throw std::runtime_error(
        "RootReader::writeData() | Can't open data file: " + OutputFilePath);

  // TTree branch variables
  TClonesArray *fParticles;
  double feventWeight;
  double feventEff;
  int fCharge;
  int fFlavour;

  auto const &Events(Data->events());

  TTree Tree(TreeName.c_str(), TreeName.c_str());
  unsigned int numPart = Events.at(0).numParticles();
  fParticles = new TClonesArray("TParticle", numPart);
  Tree.Branch("Particles", &fParticles);
  Tree.Branch("weight", &feventWeight, "weight/D");
  Tree.Branch("eff", &feventEff, "weight/D");
  Tree.Branch("charge", &fCharge, "charge/I");
  Tree.Branch("flavour", &fFlavour, "flavour/I");
  TClonesArray &partArray = *fParticles;

  for (auto const &evt : Events) {
    fParticles->Clear();
    feventWeight = evt.weight();
    fCharge = evt.charge();
    feventEff = evt.efficiency();

    TLorentzVector motherMomentum(0, 0, 0, evt.cmsEnergy());
    for (unsigned int i = 0; i < numPart; i++) {
      Particle oldParticle = evt.particle(i);
      TLorentzVector oldMomentum(oldParticle.px(), oldParticle.py(),
                                 oldParticle.pz(), oldParticle.e());
      new (partArray[i]) TParticle(oldParticle.pid(), 1, 0, 0, 0, 0,
                                   oldMomentum, motherMomentum);
    }
    Tree.Fill();
  }
  Tree.Write("", TObject::kOverwrite, 0);
  File.Close();
}

} // namespace Data
} // namespace ComPWA
