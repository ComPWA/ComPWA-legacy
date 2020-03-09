// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Data_RootDataIOTest

#include "Data/Root/RootDataIO.hpp"
#include "Core/Logging.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"

#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TTree.h"

#include <boost/test/unit_test.hpp>

#include <memory>
#include <string>

namespace ComPWA {
namespace Data {

/// Example macro to demonstrate how to generate a `TTree` that can be handled
/// by the functions in the `ComPWA::Data::Root` namespace.
void createRootTree(const std::string &OutputFileName,
                    const std::string &TreeName, int NumberOfEvents = 100) {

  TFile File(OutputFileName.c_str(), "RECREATE");
  TTree Tree(TreeName.c_str(), TreeName.c_str());
  double Weight;
  TLorentzVector *Dm = nullptr;
  TLorentzVector *Pip1 = nullptr;
  TLorentzVector *Pip2 = nullptr;
  Tree.Branch("weight", &Weight);
  Tree.Branch("_411", &Dm);
  Tree.Branch("211_1", &Pip1);
  Tree.Branch("211_2", &Pip2);

  TLorentzVector p4_cms(0, 0, 0, 3.77);
  double masses[3] = {0.5, 0.5, 0.5};
  TGenPhaseSpace gen;
  gen.SetDecay(p4_cms, 3, masses);

  LOG(INFO) << "MOCK writing vector of " << NumberOfEvents << " to file "
            << OutputFileName;
  TRandom3 rand(0);
  for (int i = 0; i < NumberOfEvents; i++) {
    Weight = rand.Uniform(0, 1);
    gen.Generate();
    Dm = gen.GetDecay(0);
    Pip1 = gen.GetDecay(1);
    Pip2 = gen.GetDecay(2);
    Tree.Fill();
  }
  Tree.Write();
  File.Close();
}

EventCollection generateSample(double InitialStateMass, std::vector<pid> Pids,
                               std::vector<double> FinalState,
                               unsigned int NumberOfEvents = 100) {
  using namespace ComPWA::Data::Root;
  RootGenerator Generator({0., 0., 0., InitialStateMass}, FinalState, Pids);
  RootUniformRealGenerator RandomGenerator(305896);
  return Generator.generate(NumberOfEvents, RandomGenerator);
}

BOOST_AUTO_TEST_SUITE(RootData);

BOOST_AUTO_TEST_CASE(SimpleWriteCheck) {
  ComPWA::Logging log("TRACE");

  const char *FileName = "RootReaderTest-WriteCheck.root";
  const char *TreeName = "tree";
  unsigned int NumberOfEvents = 100;

  using namespace ComPWA::Data::Root;
  std::vector<pid> Pids{-411, 211, 211};
  auto SampleOut = generateSample(1.864, Pids, {0.5, 0.5, 0.5}, NumberOfEvents);
  writeData(SampleOut, FileName, TreeName);

  TFile File(FileName);
  if (File.IsZombie())
    throw std::runtime_error("Could not load file");
  auto Tree = dynamic_cast<TTree *>(File.Get(TreeName));
  if (!Tree)
    throw std::runtime_error("TFile not contain expected TTree");
  if (!Tree->GetBranch("weights"))
    throw std::runtime_error("TTree not contain weights");
  if (Tree->GetListOfBranches()->GetEntries() != (Int_t)Pids.size() + 1)
    throw std::runtime_error("TTree does not contain " +
                             std::to_string(Pids.size()) + " PIDs");

  double Weight;
  std::vector<TLorentzVector *> FourMomenta(Pids.size());
  if (Tree->SetBranchAddress("weights", &Weight))
    throw std::runtime_error(
        "Could not set branch address of branch \"weights\"");
  size_t i = 0;
  TIter Next(Tree->GetListOfBranches());
  while (TObject *obj = Next()) {
    std::string BranchName = obj->GetName();
    if (BranchName == "weights")
      continue;
    if (Tree->SetBranchAddress(BranchName.c_str(), &FourMomenta.at(i)))
      throw std::runtime_error("Could not set branch address of branch \"" +
                               BranchName + "\"");
    ++i;
  }

  BOOST_CHECK_EQUAL(Tree->GetEntries(), NumberOfEvents);
  for (Long64_t i = 0; i < NumberOfEvents; ++i) {
    Tree->GetEntry(i);
    const auto &Event = SampleOut.Events.at(i);
    BOOST_CHECK_EQUAL(Weight, Event.Weight);
    for (size_t j = 0; j < Pids.size(); ++j) {
      BOOST_CHECK_EQUAL(FourMomenta.at(j)->E(), Event.FourMomenta.at(j).e());
      BOOST_CHECK_EQUAL(FourMomenta.at(j)->Px(), Event.FourMomenta.at(j).px());
      BOOST_CHECK_EQUAL(FourMomenta.at(j)->Py(), Event.FourMomenta.at(j).py());
      BOOST_CHECK_EQUAL(FourMomenta.at(j)->Pz(), Event.FourMomenta.at(j).pz());
    }
  }

  std::remove(FileName);
}

BOOST_AUTO_TEST_CASE(SimpleReadCheck) {
  ComPWA::Logging log("TRACE");

  const char *FileName = "RootReaderTest-ReadCheck.root";
  const char *TreeName = "TestTree";
  unsigned int NumberOfEvents = 10;

  using namespace ComPWA::Data::Root;
  std::vector<pid> Pids{1, 2, 3};
  auto SampleOut = generateSample(1.864, Pids, {0.5, 0.5, 0.5}, NumberOfEvents);
  writeData(SampleOut, FileName, TreeName);

  auto SampleIn = readData(FileName, TreeName);
  BOOST_CHECK_EQUAL(SampleOut.Events.size(), SampleIn.Events.size());
  for (std::size_t i = 0; i < SampleIn.Events.size(); ++i) {
    BOOST_CHECK_EQUAL(SampleIn.Events.at(i).Weight,
                      SampleOut.Events.at(i).Weight);
    BOOST_CHECK_EQUAL(SampleIn.Events.at(i).FourMomenta.front().e(),
                      SampleOut.Events.at(i).FourMomenta.front().e());
  }

  std::remove(FileName);
}

BOOST_AUTO_TEST_CASE(ReadWildcards) {
  ComPWA::Logging log("TRACE");

  using namespace ComPWA::Data::Root;
  std::vector<pid> Pids{1, 2, 3};
  auto Sample1 = generateSample(1.864, Pids, {0.5, 0.5, 0.5}, 10);
  auto Sample2 = generateSample(1.864, Pids, {0.5, 0.5, 0.5}, 20);
  auto Sample3 = generateSample(1.864, Pids, {0.5, 0.5, 0.5}, 30);
  writeData(Sample1, "RootReaderTest-writeData1.root", "Correct");
  writeData(Sample2, "RootReaderTest-writeData2.root", "WRONG");
  writeData(Sample3, "RootReaderTest-writeData3.root", "Correct");
  createRootTree("RootReaderTest-createRootTree1.root", "Correct", 15);
  createRootTree("RootReaderTest-createRootTree2.root", "WRONG", 25);
  auto SampleIn = readData("RootReaderTest-*.root", "Correct");
  BOOST_CHECK_EQUAL(SampleIn.Events.size(), 55);

  std::remove("RootReaderTest-writeData1.root");
  std::remove("RootReaderTest-writeData2.root");
  std::remove("RootReaderTest-writeData3.root");
  std::remove("RootReaderTest-createRootTree1.root");
  std::remove("RootReaderTest-createRootTree2.root");
}

BOOST_AUTO_TEST_CASE(WriteTwoTreesToSameFile) {
  ComPWA::Logging log("TRACE");

  using namespace ComPWA::Data::Root;
  std::vector<pid> Pids{1, 2, 3};
  auto SampleOut1 = generateSample(1.864, Pids, {0.5, 0.5, 0.5}, 10);
  auto SampleOut2 = generateSample(1.864, Pids, {0.5, 0.5, 0.5}, 20);

  const char *OutputFileName = "RootReaderTest-sameFile.root";
  writeData(SampleOut1, OutputFileName, "Tree1", false);
  writeData(SampleOut2, OutputFileName, "Tree2", false);

  TFile File(OutputFileName);
  if (File.IsZombie())
    throw std::runtime_error("Could not load file");
  TList *ListOfKeys = File.GetListOfKeys();
  if (!ListOfKeys)
    throw std::runtime_error("Empty list of keys");

  BOOST_CHECK_EQUAL(ListOfKeys->GetEntries(), 2);
  BOOST_CHECK_EQUAL(ListOfKeys->At(0)->GetName(), "Tree1");
  BOOST_CHECK_EQUAL(ListOfKeys->At(1)->GetName(), "Tree2");

  std::remove(OutputFileName);
}

BOOST_AUTO_TEST_CASE(OpenWrongBranches) {
  ComPWA::Logging log("TRACE");

  // Create ROOT file with wrong branch names
  const char *FileName = "WrongFormat.root";
  const char *TreeName = "data";
  TFile File(FileName, "RECREATE");
  TTree Chain(TreeName, TreeName);
  double SomeDouble;
  Chain.Branch("SomeDouble", &SomeDouble);
  for (int i = 0; i < 5; i++)
    Chain.Fill();
  Chain.Write();
  File.Close();

  BOOST_CHECK_THROW(ComPWA::Data::Root::readData(FileName, TreeName),
                    ComPWA::CorruptFile);

  std::remove(FileName);
}

BOOST_AUTO_TEST_SUITE_END();

} // namespace Data
} // namespace ComPWA
