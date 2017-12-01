// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "TLorentzVector.h"
#include "TRandom3.h"

#include "Core/DataPoint.hpp"
#include "Core/Properties.hpp"
#include "Tools/RootGenerator.hpp"

namespace ComPWA {
namespace Tools {

RootGenerator::RootGenerator(double cmsEnergy, double m1, double m2, double m3,
                             int seed)
    : nPart(3), cmsP4(0, 0, 0, cmsEnergy) {
  gRandom = new TRandom3(0);
  if (seed != -1)
    SetSeed(seed);

  masses = new Double_t[nPart];
  masses[0] = m1;
  masses[1] = m2;
  masses[2] = m3;

  TLorentzVector W(cmsP4.GetPx(), cmsP4.GetPy(), cmsP4.GetPz(), cmsP4.GetE());
  event.SetDecay(W, nPart, masses);
}

RootGenerator::RootGenerator(std::shared_ptr<PartList> partL,
                             std::vector<pid> initialS, std::vector<pid> finalS,
                             int seed) {
  gRandom = new TRandom3(0);
  if (seed != -1)
    SetSeed(seed);
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
        "RootGenerator::RootGenerator() | More than one "
        "particle in initial State! You need to specify a cms four-momentum");

  double sqrtS = FindParticle(partL, initialS.at(0)).GetMass();
  cmsP4 = FourMomentum(0, 0, 0, sqrtS);

  masses = new Double_t[nPart];
  TLorentzVector W(0.0, 0.0, 0.0, sqrtS);    //= beam + target;
  for (unsigned int t = 0; t < nPart; t++) { // particle 0 is mother particle
    masses[t] = FindParticle(partL, finalS.at(t)).GetMass();
  }
  event.SetDecay(W, nPart, masses);
};

RootGenerator::RootGenerator(std::shared_ptr<PartList> partL,
                             std::shared_ptr<Kinematics> kin, int seed) {
  gRandom = new TRandom3(0);
  if (seed != -1)
    SetSeed(seed);
  auto finalS = kin->finalState();
  auto initialS = kin->initialState();
  nPart = finalS.size();
  if (nPart < 2)
    throw std::runtime_error(
        "RootGenerator::RootGenerator() | one particle is not enough!");
  if (nPart == 2)
    LOG(info)
        << "RootGenerator::RootGenerator() | only 2 particles in the final"
           " state! There are no degrees of freedom!";

  cmsP4 = kin->initialStateFourMomentum();
  TLorentzVector W(cmsP4.GetPx(), cmsP4.GetPy(), cmsP4.GetPz(), cmsP4.GetE());

  masses = new Double_t[nPart];
  for (unsigned int t = 0; t < nPart; t++) { // particle 0 is mother particle
    masses[t] = FindParticle(partL, finalS.at(t)).GetMass();
  }
  event.SetDecay(W, nPart, masses);
};

RootGenerator *RootGenerator::Clone() { return (new RootGenerator(*this)); }

void RootGenerator::Generate(Event &evt) {
  evt.clear();
  const double weight = event.Generate();
  for (unsigned int t = 0; t < nPart; t++) {
    TLorentzVector *p = event.GetDecay(t);
    evt.addParticle(Particle(p->X(), p->Y(), p->Z(), p->E()));
  }
  evt.setWeight(weight);

#ifndef _NDEBUG
  ComPWA::FourMomentum pFour;
  double sqrtS = cmsP4.GetInvMass();

  for (int i = 0; i < evt.numParticles(); i++)
    pFour += evt.particle(i).GetFourMomentum();
  if (pFour.GetInvMass() != cmsP4.GetInvMass()) {
    // TGenPhaseSpace calculates momenta with float precision. This can lead
    // to the case that generated events are outside the available
    // phase space region. Haven't found a solution yet.
    // You can increase the numerical precison in the following compare
    // function.
    if (!ComPWA::equal(pFour.GetInvMass(), sqrtS, 100)) {
      LOG(error) << pFour.GetInvMass() << " - " << sqrtS << " = "
                 << pFour.GetInvMass() - sqrtS;
      throw std::runtime_error(
          "RootGenerator::Generate() | Invariant mass of "
          "all generate particles does not sum up to the mass of the decaying "
          "particle.");
    }
  }
#endif

  return;
}

void RootGenerator::SetSeed(unsigned int seed) { gRandom->SetSeed(seed); }

unsigned int RootGenerator::GetSeed() const { return gRandom->GetSeed(); }

double RootGenerator::GetGaussDist(double mu, double sigma) const {
  return gRandom->Gaus(mu, sigma);
}

double RootGenerator::GetUniform(double min, double max) const {
  return gRandom->Uniform(min, max);
}

void UniformTwoBodyGenerator::Generate(Event &evt) {
  double s = RootGenerator::GetUniform(minSq, maxSq);
  TLorentzVector W(0.0, 0.0, 0.0, sqrt(s)); //= beam + target;
  RootGenerator::GetGenerator()->SetDecay(W, nPart, masses);
  RootGenerator::Generate(evt);
}
} // ns::Tools
} // ns::ComPWA
