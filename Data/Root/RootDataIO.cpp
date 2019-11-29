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

class EventDataStructure {
public:
  friend std::vector<ComPWA::Event> readData(const std::string &InputFileName,
                                             const std::string &TreeName,
                                             Long64_t NumberEventsToRead);
  friend void writeData(const std::vector<ComPWA::Event> &Events,
                        const std::string &OutputFilePath,
                        const std::string &TreeName, bool OverwriteFile);
  EventDataStructure(Int_t NumberOfParticles = 0);

private:
  TClonesArray Particles;
  TClonesArray *ParticlesPointer;
  double Weight;
};

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
  EventDataStructure TheEventStruct;
  TheChain.GetBranch("Particles")->SetAutoDelete(false);
  TheChain.SetBranchAddress("Particles", &TheEventStruct.ParticlesPointer);
  TheChain.SetBranchAddress("weight", &TheEventStruct.Weight);

  std::vector<ComPWA::Event> Events;
  Events.reserve(NumberOfEventsToRead);

  for (Long64_t i = 0; i < NumberOfEventsToRead; ++i) {
    TheChain.GetEntry(i);

    TheEventStruct.Particles.Clear();
    Event TheEvent;
    TParticle *TheParticle;

    auto NumberOfParticles = TheEventStruct.Particles.GetEntriesFast();
    for (auto part = 0; part < NumberOfParticles; part++) {
      TheParticle = nullptr;
      TheParticle = (TParticle *)TheEventStruct.Particles.At(part);
      if (!TheParticle)
        continue;
      TheEvent.ParticleList.push_back(
          Particle(TheParticle->Px(), TheParticle->Py(), TheParticle->Pz(),
                   TheParticle->Energy(), TheParticle->GetPdgCode()));
    } // particle loop
    TheEvent.Weight = TheEventStruct.Weight;

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
  TTree tree(TreeName.c_str(), TreeName.c_str());
  auto TheEventStruct = EventDataStructure(Events[0].ParticleList.size());
  tree.Branch("Particles", &TheEventStruct.ParticlesPointer);
  tree.Branch("weight", &TheEventStruct.Weight, "weight/D");

  for (auto const &evt : Events) {
    TheEventStruct.Particles.Clear();
    TheEventStruct.Weight = evt.Weight;

    TLorentzVector motherMomentum(0, 0, 0, ComPWA::calculateInvariantMass(evt));
    for (unsigned int i = 0; i < evt.ParticleList.size(); ++i) {
      const Particle &x(evt.ParticleList[i]);
      auto FourMom(x.fourMomentum());
      TLorentzVector oldMomentum(FourMom.px(), FourMom.py(), FourMom.pz(),
                                 FourMom.e());
      new (TheEventStruct.Particles[i])
          TParticle(x.pid(), 1, 0, 0, 0, 0, oldMomentum, motherMomentum);
    }
    tree.Fill();
  }
  tree.Write("", TObject::kOverwrite, 0);
  TheFile.Close();
}

} // namespace Root
} // namespace Data
} // namespace ComPWA
