/*
 * RootGenerator.cpp
 *
 *  Created on: Jun 18, 2015
 *      Author: weidenka
 */

#include <stdexcept>

#include "Core/DataPoint.hpp"
#include "Core/PhysConst.hpp"
#include "Tools/RootGenerator.hpp"

namespace ComPWA {
namespace Tools {

RootGenerator::RootGenerator(int seed) {
  gRandom = new TRandom3(0);
  if (seed != -1)
    setSeed(seed);
  Kinematics *kin = Kinematics::instance();
  auto physConst = PhysConst::Instance();
  auto finalS = kin->GetFinalState();
  auto initialS = kin->GetInitialState();
  nPart = finalS.size();
  if (nPart < 2)
    throw std::runtime_error(
        "RootGenerator::RootGenerator() | one particle is not enough!");
  if (nPart == 2)
    LOG(info)
        << "RootGenerator::RootGenerator() | only 2 particles in the final"
           " state! There are no degrees of freedom!";
  if (initialS.size() != 1)
    throw std::runtime_error(
        "RootGenerator::RootGenerator() | More than one particle in initial State!");
  
  masses = new Double_t[nPart];
  TLorentzVector W(0.0, 0.0, 0.0, physConst->FindParticle(initialS.at(0)).GetMass()); //= beam + target;
  for (unsigned int t = 0; t < nPart; t++) { // particle 0 is mother particle
    masses[t] = physConst->FindParticle(finalS.at(t)).GetMass();
  }
  event.SetDecay(W, nPart, masses);
};

RootGenerator *RootGenerator::Clone() { return (new RootGenerator(*this)); }

void RootGenerator::generate(Event &evt) {
  const double weight = event.Generate();

  Event tmp(weight, "");
  for (unsigned int t = 0; t < nPart; t++) {
    TLorentzVector *p = event.GetDecay(t);
    tmp.addParticle(Particle(p->X(), p->Y(), p->Z(), p->E()));
  }
  evt = tmp;
  return;
}

void RootGenerator::setSeed(unsigned int seed) { gRandom->SetSeed(seed); }

unsigned int RootGenerator::getSeed() { return gRandom->GetSeed(); }

double RootGenerator::getGaussDist(double mu, double sigma) {
  return gRandom->Gaus(mu, sigma);
}

double RootGenerator::getUniform() { return gRandom->Uniform(0, 1); }

void UniformTwoBodyGenerator::generate(Event &evt) {
  double s = RootGenerator::getUniform() * (maxSq - minSq) + minSq;
  TLorentzVector W(0.0, 0.0, 0.0, sqrt(s)); //= beam + target;
  RootGenerator::getGenerator()->SetDecay(W, nPart, masses);
  RootGenerator::generate(evt);
}
} /* namespace Tools */
} /* namespace ComPWA */
