// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Root-Headers
#include "RootDataIO.hpp"

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
#include "Data/DataSet.hpp"

namespace ComPWA {
namespace Data {

RootDataIO::RootDataIO(const std::string TreeName_, int NumberEventsToProcess_)
    : TreeName(TreeName_), NumberEventsToProcess(NumberEventsToProcess_) {}

std::shared_ptr<DataSet>
RootDataIO::readData(const std::string &InputFilePath) const {
  TFile File(InputFilePath.c_str());
  if (File.IsZombie())
    throw std::runtime_error("RootDataIO::RootDataIO() | "
                             "Can't open data file: " +
                             InputFilePath);

  TTree *fTree = (TTree *)File.Get(TreeName.c_str());

  if (!fTree)
    throw std::runtime_error("RootDataIO::RootDataIO() | Tree \"" + TreeName +
                             "\" can not be opened from file " + InputFilePath +
                             "! ");

  // TTree branch variables
  TClonesArray Particles("TParticle");
  TClonesArray *pParticles(&Particles);
  double feventWeight;
  int fCharge;
  int fFlavour;

  fTree->GetBranch("Particles")->SetAutoDelete(false);
  fTree->SetBranchAddress("Particles", &pParticles);
  fTree->SetBranchAddress("weight", &feventWeight);
  fTree->SetBranchAddress("charge", &fCharge);
  fTree->SetBranchAddress("flavour", &fFlavour);

  unsigned int NumberEventsToRead(NumberEventsToProcess);
  if (NumberEventsToProcess <= 0 || NumberEventsToProcess > fTree->GetEntries())
    NumberEventsToRead = fTree->GetEntries();

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
      evt.ParticleList.push_back(Particle(inN.X(), inN.Y(), inN.Z(), inN.E(),
                                          partN->GetPdgCode(), charge));
    } // particle loop
    evt.Weight = feventWeight;

    Events.push_back(evt);
  } // end event loop
  File.Close();

  return std::make_shared<DataSet>(Events);
}

void RootDataIO::writeData(std::shared_ptr<const DataSet> DataSample,
                           const std::string &OutputFilePath) const {

  LOG(INFO) << "RootDataIO::writeData(): writing current "
               "vector of events to file "
            << OutputFilePath;

  auto Events = DataSample->getEventList();

  if (0 == Events.size()) {
    LOG(ERROR) << "RootDataIO::writeData(): no events given!";
    return;
  }

  TFile File(OutputFilePath.c_str(), "RECREATE");
  if (File.IsZombie())
    throw std::runtime_error("RootDataIO::writeData(): can't open data file: " +
                             OutputFilePath);

  // TTree branch variables
  TClonesArray *fParticles;
  double feventWeight;
  int fFlavour;

  TTree Tree(TreeName.c_str(), TreeName.c_str());
  unsigned int numPart = Events[0].ParticleList.size();
  fParticles = new TClonesArray("TParticle", numPart);
  Tree.Branch("Particles", &fParticles);
  Tree.Branch("weight", &feventWeight, "weight/D");
  Tree.Branch("flavour", &fFlavour, "flavour/I");
  TClonesArray &partArray = *fParticles;

  for (auto const &evt : Events) {
    fParticles->Clear();
    feventWeight = evt.Weight;

    TLorentzVector motherMomentum(0, 0, 0, ComPWA::calculateInvariantMass(evt));
    for (unsigned int i = 0; i < evt.ParticleList.size(); ++i) {
      const Particle &x(evt.ParticleList[i]);
      TLorentzVector oldMomentum(x.px(), x.py(), x.pz(), x.e());
      new (partArray[i])
          TParticle(x.pid(), 1, 0, 0, 0, 0, oldMomentum, motherMomentum);
    }
    Tree.Fill();
  }
  Tree.Write("", TObject::kOverwrite, 0);
  File.Close();
}

} // namespace Data
} // namespace ComPWA
