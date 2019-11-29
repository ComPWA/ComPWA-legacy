// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <list>

#include "RootDataIO.hpp"

#include "TChain.h"
#include "TClonesArray.h"
#include "TError.h" // for ignoring ROOT warnings
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

std::vector<ComPWA::Event> readData(const std::string &InputFilePath,
                                    const std::string &TreeName,
                                    long long NumberOfEventsToRead) {
  // Ignore custom streamer warning and error message for missing trees
  gErrorIgnoreLevel = kFatal;

  // Use TChain to add files through wildcards if necessary
  TChain TheChain(TreeName.c_str());
  TheChain.Add(InputFilePath.c_str());

  // Test TChain quality
  auto ListOfFiles = TheChain.GetListOfFiles();
  if (!ListOfFiles || !ListOfFiles->GetEntries())
    throw std::runtime_error("Root::readData() | Unable to load files: " +
                             InputFilePath);
  if (!TheChain.GetEntries())
    throw std::runtime_error("Root::readData() | TTree \"" + TreeName +
                             "\" can not be opened from file(s) " +
                             InputFilePath + "!");
  if (NumberOfEventsToRead <= 0 || NumberOfEventsToRead > TheChain.GetEntries())
    NumberOfEventsToRead = TheChain.GetEntries();
  if (!TheChain.GetBranch("Particles") || !TheChain.GetBranch("weight"))
    throw std::runtime_error("Root::readData() | TTree does not have a "
                             "Particles and/or weight branch");

  // Set branch addresses
  TClonesArray Particles("TParticle");
  TClonesArray *ParticlesPointer(&Particles);
  double Weight;
  TheChain.GetBranch("Particles")->SetAutoDelete(false);
  TheChain.SetBranchAddress("Particles", &ParticlesPointer);
  TheChain.SetBranchAddress("weight", &Weight);

  std::vector<ComPWA::Event> Events;
  Events.reserve(NumberOfEventsToRead);

  for (Long64_t i = 0; i < NumberOfEventsToRead; ++i) {
    Particles.Clear();
    TheChain.GetEntry(i);

    Event TheEvent;

    auto NumberOfParticles = Particles.GetEntries();
    for (auto part = 0; part < NumberOfParticles; part++) {
      auto TheParticle = dynamic_cast<TParticle *>(Particles.At(part));
      if (!TheParticle)
        continue;
      TheEvent.ParticleList.push_back(
          Particle(TheParticle->Px(), TheParticle->Py(), TheParticle->Pz(),
                   TheParticle->Energy(), TheParticle->GetPdgCode()));
    } // end particle loop
    TheEvent.Weight = Weight;

    Events.push_back(TheEvent);
  } // end event loop

  return Events;
}

void writeData(const std::vector<ComPWA::Event> &Events,
               const std::string &OutputFileName, const std::string &TreeName,
               bool OverwriteFile) {
  // Ignore custom streamer warning
  gErrorIgnoreLevel = kError;

  if (0 == Events.size()) {
    LOG(ERROR) << "Root::writeData(): no events given!";
    return;
  }

  LOG(INFO) << "Root::writeData(): writing vector of " << Events.size()
            << " events to file " << OutputFileName;

  std::string WriteFlag{"UPDATE"};
  if (OverwriteFile)
    WriteFlag = "RECREATE";

  TFile TheFile(OutputFileName.c_str(), "UPDATE");
  if (TheFile.IsZombie())
    throw std::runtime_error("Root::writeData(): can't open data file: " +
                             OutputFileName);

  // TTree branch variables
  TTree TheTree(TreeName.c_str(), TreeName.c_str());
  auto ParticlesPointer =
      new TClonesArray("TParticle", Events[0].ParticleList.size());
  double Weight;
  TheTree.Branch("Particles", &ParticlesPointer);
  TheTree.Branch("weight", &Weight, "weight/D");
  auto &Particles = *ParticlesPointer;

  for (auto const &Event : Events) {
    Particles.Clear();
    Weight = Event.Weight;
    double Mass = ComPWA::calculateInvariantMass(Event);
    TLorentzVector motherMomentum(0, 0, 0, Mass);
    for (unsigned int i = 0; i < Event.ParticleList.size(); ++i) {
      auto &Particle = Event.ParticleList[i];
      auto FourMom(Particle.fourMomentum());
      TLorentzVector oldMomentum(FourMom.px(), FourMom.py(), FourMom.pz(),
                                 FourMom.e());
      Particles[i] = new TParticle(Particle.pid(), 1, 0, 0, 0, 0, oldMomentum,
                                   motherMomentum);
      // TClonesArray owns its objects
    }
    TheTree.Fill();
  }
  TheTree.Write("", TObject::kOverwrite, 0);
  TheFile.Close();
}

} // namespace Root
} // namespace Data
} // namespace ComPWA
