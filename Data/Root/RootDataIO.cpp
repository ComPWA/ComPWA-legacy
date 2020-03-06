// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "RootDataIO.hpp"

#include "TChain.h"
#include "TClonesArray.h"
#include "TError.h" // for ignoring ROOT warnings
#include "TFile.h"
#include "TLorentzVector.h"
#include "TParticle.h"

#include "Core/Exceptions.hpp"
#include "Core/Generator.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"

namespace ComPWA {
namespace Data {
namespace Root {

std::vector<std::string> pidsToUniqueStrings(std::vector<pid> Pids) {
  std::vector<std::string> PidStrings;
  for (auto Pid : Pids) {
    auto PidString = std::to_string(Pid);
    if (PidString.front() == '-')
      PidString.front() = '_'; // needed because TBrowser cannot handle minus
    PidStrings.push_back(PidString);
  }
  for (const auto PidString : PidStrings) {
    int Count = 1;
    if (std::count(PidStrings.begin(), PidStrings.end(), PidString) > 1) {
      auto Iter = PidStrings.begin();
      while (Iter != PidStrings.end()) {
        Iter = std::find(Iter, PidStrings.end(), PidString);
        *Iter += "_" + std::to_string(Count);
        ++Count;
      }
    }
  }
  return PidStrings;
}

std::vector<pid> uniqueStringsToPids(std::vector<std::string> UniqueStrings) {
  std::vector<pid> Pids;
  for (const auto &PidString : UniqueStrings) {
    try {
      auto StrippedPidString = PidString.substr(0, PidString.find("_", 0));
      Pids.push_back(std::stoi(StrippedPidString));
    } catch (const std::invalid_argument &e) {
      throw ComPWA::CorruptFile("Branches \"" + PidString +
                                "\" cannot be converted to a PID");
    }
  }
  return Pids;
}

ComPWA::EventCollection readData(const std::string &InputFilePath,
                                 const std::string &TreeName,
                                 long long NumberOfEventsToRead) {
  /// -# Ignore custom streamer warning and error message for missing trees
  auto temp_ErrorIgnoreLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kBreak;

  /// -# Use TChain to add files through wildcards if necessary
  TChain Chain(TreeName.c_str());
  Chain.Add(InputFilePath.c_str());

  /// -# Test TChain quality
  auto ListOfFiles = Chain.GetListOfFiles();
  if (!ListOfFiles || !ListOfFiles->GetEntries()) {
    throw ComPWA::BadConfig("Root::readData() | Unable to load files: " +
                            InputFilePath);
  }
  if (!Chain.GetEntries()) {
    throw ComPWA::CorruptFile("Root::readData() | TTree \"" + TreeName +
                              "\" cannot be opened from file(s) " +
                              InputFilePath + "!");
  }
  if (NumberOfEventsToRead <= 0 || NumberOfEventsToRead > Chain.GetEntries()) {
    NumberOfEventsToRead = Chain.GetEntries();
  }
  if (!Chain.GetBranch("weights")) {
    throw ComPWA::CorruptFile("Root::readData() | TTree \"" + TreeName +
                              "\" in file \"" + InputFilePath +
                              "\" does not contain weights branch");
  }

  /// -# Get PIDs
  std::vector<std::string> PidStrings;
  auto Branches = Chain.GetListOfBranches();
  for (int i = 0; i < Branches->GetEntries(); ++i) {
    auto Branch = Branches->At(i);
    if (!Branch)
      continue;
    std::string Name = Branch->GetName();
    if (Name == "weights")
      continue;
    PidStrings.push_back(Name);
  }
  if (!PidStrings.size()) {
    throw ComPWA::CorruptFile("Root::readData() | TTree \"" + TreeName +
                              "\" in file \"" + InputFilePath +
                              "\" does not contain any LorentzVector branches");
  }
  EventCollection ImportedDataSample{uniqueStringsToPids(PidStrings)};

  /// -# Set branch addresses
  double Weight;
  std::vector<TLorentzVector *> LorentzVectors(ImportedDataSample.Pids.size());
  if (Chain.SetBranchAddress("weights", &Weight))
    throw std::runtime_error(
        "Could not set branch address of branch \"weights\"");
  size_t i = 0;
  TIter Next(Chain.GetListOfBranches());
  while (TObject *obj = Next()) {
    std::string BranchName = obj->GetName();
    if (BranchName == "weights")
      continue;
    if (Chain.SetBranchAddress(BranchName.c_str(), &LorentzVectors.at(i)))
      throw std::runtime_error("Could not set branch address of branch \"" +
                               BranchName + "\"");
    ++i;
  }

  /// -# Import data sample
  ImportedDataSample.Events.resize(NumberOfEventsToRead);
  for (Long64_t i = 0; i < NumberOfEventsToRead; ++i) {
    Chain.GetEntry(i);
    auto &Event = ImportedDataSample.Events.at(i);
    Event.Weight = Weight;
    for (const auto &LorentzVector : LorentzVectors) {
      Event.FourMomenta.push_back(
          FourMomentum(LorentzVector->Px(), LorentzVector->Py(),
                       LorentzVector->Pz(), LorentzVector->Energy()));
    }
  }

  gErrorIgnoreLevel = temp_ErrorIgnoreLevel;
  return ImportedDataSample;
}

void writeData(const EventCollection &OutputSample,
               const std::string &OutputFileName, const std::string &TreeName,
               bool OverwriteFile) {
  /// -# Ignore custom streamer warning
  auto temp_ErrorIgnoreLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kBreak;

  /// -# Check EventCollection quality
  if (OutputSample.Events.size() == 0) {
    throw ComPWA::CorruptFile("Root::writeData(): no events given!");
  }

  if (!OutputSample.checkPidMatchesEvents()) {
    throw ComPWA::CorruptFile(
        "Root::writeData(): number of PIDs in EventCollection "
        "does not match that of the PID list!");
  }

  /// -# Open or create ROOT file
  std::string WriteFlag{"UPDATE"};
  if (OverwriteFile)
    WriteFlag = "RECREATE";
  TFile File(OutputFileName.c_str(), WriteFlag.c_str());
  if (File.IsZombie()) {
    throw std::runtime_error("Root::writeData(): can't open data file: " +
                             OutputFileName);
  }
  LOG(INFO) << "Root::writeData(): writing vector of "
            << OutputSample.Events.size() << " events to file "
            << OutputFileName;

  /// -# Define TTree its branches
  TTree Tree(TreeName.c_str(), TreeName.c_str());
  auto BranchNames = pidsToUniqueStrings(OutputSample.Pids);
  std::vector<TLorentzVector> LorentzVectors(OutputSample.Pids.size());
  double Weight;
  Tree.Branch("weights", &Weight);
  for (size_t i = 0; i < LorentzVectors.size(); ++i) {
    Tree.Branch(BranchNames.at(i).c_str(), &LorentzVectors.at(i));
  }
  /// -# Fill tree
  for (auto const &Event : OutputSample.Events) {
    Weight = Event.Weight;
    for (size_t i = 0; i < Event.FourMomenta.size(); ++i) {
      auto &LorentzVector = LorentzVectors.at(i);
      LorentzVector.SetPx(Event.FourMomenta.at(i).px());
      LorentzVector.SetPy(Event.FourMomenta.at(i).py());
      LorentzVector.SetPz(Event.FourMomenta.at(i).pz());
      LorentzVector.SetE(Event.FourMomenta.at(i).e());
    }
    Tree.Fill();
  }
  Tree.Write();
  File.Close();
  gErrorIgnoreLevel = temp_ErrorIgnoreLevel;
}

} // namespace Root
} // namespace Data
} // namespace ComPWA
