// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Data_RootDataIOTest

#include "Data/Root/RootDataIO.hpp"
#include "Core/Logging.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"

#include "TClonesArray.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TParticle.h"
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

  TFile TheFile(OutputFileName.c_str(), "RECREATE");
  TTree TheTree(TreeName.c_str(), TreeName.c_str());
  TClonesArray parts("TParticle", 3);
  double weight;
  TheTree.Branch("Particles", &parts);
  TheTree.Branch("weight", &weight);

  TLorentzVector p4_cms(0, 0, 0, 3.77);
  double masses[3] = {0.5, 0.5, 0.5};
  TGenPhaseSpace gen;
  gen.SetDecay(p4_cms, 3, masses);

  TRandom3 rand(0);
  for (int i = 0; i < NumberOfEvents; i++) {
    weight = rand.Uniform(0, 1);
    parts = TClonesArray("TParticle", 3);
    gen.Generate();

    TLorentzVector p4_1 = *gen.GetDecay(0);
    TLorentzVector p4_2 = *gen.GetDecay(1);
    TLorentzVector p4_3 = *gen.GetDecay(2);
    new (parts[0]) TParticle(310, 1, 0, 0, 0, 0, p4_1, p4_cms);
    new (parts[1]) TParticle(-321, 1, 0, 0, 0, 0, p4_2, p4_cms);
    new (parts[2]) TParticle(321, 1, 0, 0, 0, 0, p4_3, p4_cms);
    TheTree.Fill();
  }
  TheTree.Write();
  TheFile.Close();
}

std::vector<ComPWA::Event> generateSample(double InitialStateMass,
                                          std::vector<double> FinalState,
                                          unsigned int NumberOfEvents = 100) {
  using namespace ComPWA::Data::Root;
  RootGenerator gen({0., 0., 0., InitialStateMass}, FinalState);
  RootUniformRealGenerator RandomGenerator(305896);
  return ComPWA::Data::generatePhsp(NumberOfEvents, gen, RandomGenerator);
}

BOOST_AUTO_TEST_SUITE(RootData);

BOOST_AUTO_TEST_CASE(SimpleWriteCheck) {
  ComPWA::Logging log("TRACE");

  const char *FileName = "RootReaderTest-WriteCheck.root";
  const char *TreeName = "tree";
  unsigned int NumberOfEvents = 10;

  using namespace ComPWA::Data::Root;
  auto SampleOut = generateSample(1.864, {0.5, 0.5, 0.5}, NumberOfEvents);
  writeData(SampleOut, FileName, TreeName);

  TFile TheFile(FileName);
  if (TheFile.IsZombie())
    throw std::runtime_error("Could not load file");
  auto TheTree = dynamic_cast<TTree *>(TheFile.Get(TreeName));
  if (!TheTree)
    throw std::runtime_error("TFile not contain expected TTree");
  if (!TheTree->GetBranch("Particles") || !TheTree->GetBranch("weight"))
    throw std::runtime_error("TTree not contain expected branches");

  TClonesArray Particles("TParticle");
  auto ParticlesPointer(&Particles); // ROOT business...
  double Weight;
  TheTree->SetBranchAddress("Particles", &ParticlesPointer);
  TheTree->SetBranchAddress("weight", &Weight);

  BOOST_CHECK_EQUAL(TheTree->GetEntries(), NumberOfEvents);
  for (unsigned int i = 0; i < NumberOfEvents; ++i) {
    TheTree->GetEntry(i);
    BOOST_CHECK_EQUAL(Weight, SampleOut[i].Weight);
    BOOST_CHECK_EQUAL(Particles.GetEntries(), SampleOut[i].ParticleList.size());
  }

  std::remove(FileName);
}

BOOST_AUTO_TEST_CASE(SimpleReadCheck) {
  ComPWA::Logging log("TRACE");

  const char *FileName = "RootReaderTest-ReadCheck.root";
  const char *TreeName = "TestTree";
  unsigned int NumberOfEvents = 10;

  using namespace ComPWA::Data::Root;
  auto SampleOut = generateSample(1.864, {0.5, 0.5, 0.5}, NumberOfEvents);
  writeData(SampleOut, FileName, TreeName);

  auto SampleIn(readData(FileName, TreeName));
  BOOST_CHECK_EQUAL(SampleOut.size(), SampleIn.size());
  for (std::size_t i = 0; i < SampleIn.size(); ++i) {
    BOOST_CHECK_EQUAL(SampleIn[i].Weight, SampleOut[i].Weight);
    BOOST_CHECK_EQUAL(SampleIn[i].ParticleList.size(),
                      SampleOut[i].ParticleList.size());
    BOOST_CHECK_EQUAL(SampleIn[i].ParticleList.front().fourMomentum().e(),
                      SampleOut[i].ParticleList.front().fourMomentum().e());
  }

  std::remove(FileName);
}

BOOST_AUTO_TEST_CASE(ReadWildcards) {
  ComPWA::Logging log("TRACE");

  using namespace ComPWA::Data::Root;
  auto Sample1 = generateSample(1.864, {0.5, 0.5, 0.5}, 10);
  auto Sample2 = generateSample(1.864, {0.5, 0.5, 0.5}, 20);
  auto Sample3 = generateSample(1.864, {0.5, 0.5, 0.5}, 30);
  writeData(Sample1, "RootReaderTest-writeData1.root", "Correct");
  writeData(Sample2, "RootReaderTest-writeData2.root", "WRONG");
  writeData(Sample3, "RootReaderTest-writeData3.root", "Correct");
  createRootTree("RootReaderTest-createRootTree1.root", "Correct", 15);
  createRootTree("RootReaderTest-createRootTree2.root", "WRONG", 25);

  auto SampleIn(readData("RootReaderTest-*.root", "Correct"));
  BOOST_CHECK_EQUAL(SampleIn.size(), 55);

  std::remove("RootReaderTest-writeData1.root");
  std::remove("RootReaderTest-writeData2.root");
  std::remove("RootReaderTest-writeData3.root");
  std::remove("RootReaderTest-createRootTree1.root");
  std::remove("RootReaderTest-createRootTree2.root");
}

BOOST_AUTO_TEST_CASE(WriteTwoTreesToSameFile) {
  ComPWA::Logging log("TRACE");

  using namespace ComPWA::Data::Root;
  auto SampleOut1 = generateSample(1.864, {0.5, 0.5, 0.5}, 10);
  auto SampleOut2 = generateSample(1.864, {0.5, 0.5, 0.5}, 20);

  const char *OutputFileName = "RootReaderTest-sameFile.root";
  writeData(SampleOut1, OutputFileName, "Tree1", false);
  writeData(SampleOut2, OutputFileName, "Tree2", false);

  TFile TheFile(OutputFileName);
  if (TheFile.IsZombie())
    throw std::runtime_error("Could not load file");
  TList *ListOfKeys = TheFile.GetListOfKeys();
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
  TFile TheFile(FileName, "RECREATE");
  TTree TheChain(TreeName, TreeName);
  double SomeDouble;
  TheChain.Branch("SomeDouble", &SomeDouble);
  for (int i = 0; i < 5; i++)
    TheChain.Fill();
  TheChain.Write();
  TheFile.Close();

  try {
    ComPWA::Data::Root::readData(FileName, TreeName);
  } catch (const std::runtime_error &e) {
    BOOST_CHECK_EQUAL(e.what(), "Root::readData() | TTree does not have a "
                                "Particles and/or weight branch");
  }

  std::remove(FileName);
}

BOOST_AUTO_TEST_SUITE_END();

} // namespace Data
} // namespace ComPWA
