//-------------------------------------------------------------------------------
// Copyright (c) 2016 Mathias Michel.
//
// This file is part of ComPWA
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <utility>
#include <stdexcept>

#include "DataReader/JakeReader/JakeReader.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Generator.hpp"

#include "TParticle.h"

#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
using namespace boost::log;

namespace ComPWA {
namespace DataReader {
namespace JakeReader {

JakeReader::JakeReader() {
  fFile = 0; // need to do this to avoid seg. violation when destructor is
             // called
}

JakeReader::JakeReader(TTree *tr, int size, const bool binned)
    : readSize(size), fBinned(binned) {
  fTree = tr;
  fFile = 0; // need to do this to avoid seg. violation when destructor is
             // called
  read();
}

JakeReader::JakeReader(const std::string inRootFile,
                       const std::string inTreeName, int size,
                       const bool binned)
    : readSize(size), fileName(inRootFile), treeName(inTreeName),
      fBinned(binned) {
  fFile = new TFile(fileName.c_str());
  if (fFile->IsZombie())
    throw std::runtime_error("RootReader::RootReader() can't open data file: " +
                             inRootFile);
  fTree = (TTree *)fFile->Get(treeName.c_str());
  if (!fTree)
    throw std::runtime_error("RootReader::RootReader() Tree \"" + treeName +
                             "\" can not be opened can't from file " +
                             inRootFile + "! ");
  read();
  fFile->Close();
}

JakeReader::~JakeReader() {}

bool JakeReader::hasWeights() {
  bool has = 0;
  for (unsigned int evt = 0; evt < fEvents.size(); evt++) {
    if (fEvents.at(evt).getWeight() != 1.) {
      has = 1;
      break;
    }
  }
  return has;
}

JakeReader *JakeReader::Clone() const { return new JakeReader(*this); }

JakeReader *JakeReader::EmptyClone() const { return new JakeReader(); }

void JakeReader::Clear() { fEvents.clear(); }
void JakeReader::setEfficiency(std::shared_ptr<Efficiency> eff) {
  for (unsigned int evt = 0; evt < fEvents.size(); evt++) {
    dataPoint e(fEvents.at(evt));
    double val = eff->evaluate(e);
    fEvents.at(evt).setEfficiency(val);
  }
}
void JakeReader::resetEfficiency(double e) {
  for (unsigned int evt = 0; evt < fEvents.size(); evt++) {
    fEvents.at(evt).setEfficiency(e);
  }
}

void JakeReader::read() {
  // fParticles = new TClonesArray("TParticle");
  // fTree->GetBranch("kin")->SetAutoDelete(false);
  fTree->SetBranchAddress("e", fe);
  fTree->SetBranchAddress("px", fpx);
  fTree->SetBranchAddress("py", fpy);
  fTree->SetBranchAddress("pz", fpz);
  fTree->SetBranchAddress("weight", &feventWeight);
  //	fEvent=0;
  bin();
  storeEvents();
}

const std::vector<std::string> &JakeReader::getVariableNames() {
  if (!fVarNames.size()) { // TODO: init
    fVarNames.push_back("dataname1");
    fVarNames.push_back("dataname2");
  }
  return fVarNames;
}
void JakeReader::resetWeights(double w) {
  for (unsigned int i = 0; i < fEvents.size(); i++)
    fEvents.at(i).setWeight(w);
  return;
}

Event &JakeReader::getEvent(const int i) { return fEvents.at(i); }

const int JakeReader::getBin(const int i, double &m12, double &weight) {
  if (!fBinned)
    return -1;

  m12 = fBins[i].first;
  weight = fBins[i].second;

  return 1;
}

void JakeReader::storeEvents() {
  if (readSize <= 0 || readSize > fTree->GetEntries())
    readSize = fTree->GetEntries();
  for (unsigned int evt = 0; evt < readSize; evt++) {
    Event tmp;
    // fParticles->Clear();
    fTree->GetEntry(evt);

    // Get number of particle in TClonesrray
    unsigned int nParts = 3;

    // TParticle* partN;
    // TLorentzVector inN;
    for (unsigned int part = 0; part < nParts; part++) {
      // partN = 0;
      // partN = (TParticle*) fParticles->At(part);
      // if(!partN) continue;
      // partN->Momentum(inN);
      tmp.addParticle(Particle(fpx[part], fpy[part], fpz[part], fe[part], 0));
    } // particle loop
    tmp.setWeight(feventWeight);
    tmp.setCharge(0);
    tmp.setFlavour(0);
    tmp.setEfficiency(1.);

    fEvents.push_back(tmp);
  } // event loop
}

void JakeReader::writeData(std::string file, std::string trName) {
  if (file != "")
    fileName = file;
  if (trName != "")
    treeName = trName;
  /* LOG(info) << "RootReader: writing current vector of events to file
   "<<fileName;
   TFile* fFile = new TFile(fileName.c_str(),"UPDATE");
   if(fFile->IsZombie())
          throw std::runtime_error("RootReader::RootReader() can't open data
   file: "+fileName);

   fTree = new TTree(treeName.c_str(),treeName.c_str());
   unsigned int numPart = fEvents[0].getNParticles();
   fParticles = new TClonesArray("TParticle",numPart);
   fTree->Branch("Particles",&fParticles);
   fTree->Branch("weight",&feventWeight,"weight/D");
   fTree->Branch("eff",&feventEff,"weight/D");
   fTree->Branch("charge",&fCharge,"charge/I");
   fTree->Branch("flavour",&fFlavour,"flavour/I");
   TClonesArray &partArray = *fParticles;

   TLorentzVector motherMomentum(0,0,0,Kinematics::instance()->getMotherMass());
   for(std::vector<Event>::iterator it=fEvents.begin(); it!=fEvents.end();
   it++){
          fParticles->Clear();
          feventWeight = (*it).getWeight();
          fCharge= (*it).getCharge();
          fFlavour= (*it).getFlavour();
          feventEff = (*it).getEfficiency();
          for(unsigned int i=0; i<numPart; i++){
                 const Particle oldParticle = (*it).getParticle(i);
                 TLorentzVector
   oldMomentum(oldParticle.px,oldParticle.py,oldParticle.pz,oldParticle.E);
                 new(partArray[i])
   TParticle(oldParticle.pid,1,0,0,0,0,oldMomentum,motherMomentum);
          }
          fTree->Fill();
   }
   fTree->Write("",TObject::kOverwrite,0);
   fFile->Close();*/
  return;
}

void JakeReader::bin() {
  fmaxBins = 500; // TODO setter, consturctor
  TLorentzVector in1, in2;

  // initialize min & max
  fTree->GetEntry(0);
  unsigned int nParts = 3;
  if (nParts != 2)
    return;
  /*TParticle* part1 = (TParticle*) fParticles->At(0); //pip
  TParticle* part2 = (TParticle*) fParticles->At(1); //pim
  if(!part1 || !part2) return;
  part1->Momentum(in1);
  part2->Momentum(in2);
  double inm12=(in1+in2).Mag2();
  min = max = inm12;

  //find min and max
  for(unsigned int evt=1; evt<fTree->GetEntries(); evt++){
         fTree->GetEntry(evt);
         unsigned int nParts = fParticles->GetEntriesFast();
         if(nParts!=2) return;

         part1 = (TParticle*) fParticles->At(0); //pip
         part2 = (TParticle*) fParticles->At(1); //pim
         if(!part1 || !part2) return;
         part1->Momentum(in1);
         part2->Momentum(in2);

         inm12=(in1+in2).Mag2();

         if( min>inm12 ) min=inm12;
         if( max<inm12 ) max=inm12;
  }

  //initialize bins with weight zero
  min = sqrt(min);
  max = sqrt(max);
  double step = (max-min)/(double)fmaxBins;
  for(unsigned int bin=0; bin<fmaxBins; bin++){
         fBins[bin] = std::make_pair(min+bin*step,0.0);
  }

  //fill bins
  for(unsigned int evt=0; evt<fmaxBins; evt++){
         fTree->GetEntry(evt);
         unsigned int nParts = fParticles->GetEntriesFast();
               if(nParts!=2) return;

               TParticle* part1 = (TParticle*) fParticles->At(0); //pip
               TParticle* part2 = (TParticle*) fParticles->At(1); //pim
               if(!part1 || !part2) return;
               part1->Momentum(in1);
               part2->Momentum(in2);
               inm12=sqrt((in1+in2).Mag2());

               int bin = (int)((inm12-min)/step);
               fBins[bin].second += 1.;
       }*/
}

std::vector<dataPoint> JakeReader::getDataPoints() {
  std::vector<dataPoint> vecPoint;
  for (int i = 0; i < fEvents.size(); i++)
    vecPoint.push_back(dataPoint(fEvents.at(i)));
  return vecPoint;
}

} /* namespace JakeReader */
} /* namespace DataReader */
} /* namespace ComPWA */
