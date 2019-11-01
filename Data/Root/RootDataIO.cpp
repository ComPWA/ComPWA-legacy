// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Root-Headers
#include "RootDataIO.hpp"

#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include "Core/Generator.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"

namespace ComPWA {
namespace Data {
namespace Root {

RootDataIO::RootDataIO(const std::string TreeName_,
                       std::size_t NumberEventsToProcess_)
    : TreeName(TreeName_), NumberEventsToProcess(NumberEventsToProcess_) {}

std::vector<ComPWA::Event>
RootDataIO::readData(const std::string &InputFilePath) const {
  /// @todo Use of `TChain` in this way may result in run-time errors for larger
  /// data samples.
  TChain chain(TreeName.c_str());
  chain.Add(InputFilePath.c_str());
  if (!chain.GetListOfFiles()->GetEntriesFast())
    throw std::runtime_error("RootDataIO::RootDataIO() | "
                             "Unable to load files: " +
                             InputFilePath);

  if (!chain.GetEntriesFast())
    throw std::runtime_error("RootDataIO::RootDataIO() | Tree \"" + TreeName +
                             "\" can not be opened from file " + InputFilePath +
                             "!");

  // TTree branch variables
  TClonesArray Particles("TParticle");
  TClonesArray *pParticles(&Particles);
  double feventWeight;

  chain.GetBranch("Particles")->SetAutoDelete(false);
  chain.SetBranchAddress("Particles", &pParticles);
  chain.SetBranchAddress("weight", &feventWeight);

  unsigned int NumberEventsToRead(NumberEventsToProcess);
  if (NumberEventsToProcess <= 0 || NumberEventsToProcess > chain.GetEntries())
    NumberEventsToRead = chain.GetEntries();

  std::vector<ComPWA::Event> Events;
  Events.reserve(NumberEventsToRead);

  for (unsigned int i = 0; i < NumberEventsToRead; ++i) {
    Event evt;
    Particles.Clear();
    chain.GetEntry(i);

    // Get number of particle in TClonesrray
    auto nParts = Particles.GetEntriesFast();

    TParticle *partN;
    TLorentzVector inN;
    for (auto part = 0; part < nParts; part++) {
      partN = 0;
      partN = (TParticle *)Particles.At(part);
      if (!partN)
        continue;
      partN->Momentum(inN);
      evt.ParticleList.push_back(
          Particle(inN.X(), inN.Y(), inN.Z(), inN.E(), partN->GetPdgCode()));
    } // particle loop
    evt.Weight = feventWeight;

    Events.push_back(evt);
  } // end event loop

  return Events;
}

void RootDataIO::writeData(const std::vector<ComPWA::Event> &Events,
                           const std::string &OutputFilePath) const {

  LOG(INFO) << "RootDataIO::writeData(): writing current "
               "vector of events to file "
            << OutputFilePath;

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
      auto FourMom(x.fourMomentum());
      TLorentzVector oldMomentum(FourMom.px(), FourMom.py(), FourMom.pz(),
                                 FourMom.e());
      new (partArray[i])
          TParticle(x.pid(), 1, 0, 0, 0, 0, oldMomentum, motherMomentum);
    }
    Tree.Fill();
  }
  Tree.Write("", TObject::kOverwrite, 0);
  File.Close();
}

} // namespace Root
} // namespace Data
} // namespace ComPWA
