

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

RootGenerator::RootGenerator(double cmsEnergy, double m1, double m2, double m3,
                             int seed) {
  gRandom = new TRandom3(0);
  if (seed != -1)
    SetSeed(seed);

  nPart = 3;

  masses = new Double_t[nPart];
  masses[0] = m1;
  masses[1] = m2;
  masses[2] = m3;

  sqrtS = cmsEnergy;
  TLorentzVector W(0.0, 0.0, 0.0, sqrtS);
  event.SetDecay(W, nPart, masses);
}

RootGenerator::RootGenerator(std::vector<pid> initialS, std::vector<pid> finalS,
                             int seed) {
  gRandom = new TRandom3(0);
  if (seed != -1)
    SetSeed(seed);
  auto physConst = PhysConst::Instance();
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
        "particle in initial State! Currently we can not deal with that!");

  sqrtS = physConst->FindParticle(initialS.at(0)).GetMass();

  masses = new Double_t[nPart];
  TLorentzVector W(0.0, 0.0, 0.0, sqrtS);    //= beam + target;
  for (unsigned int t = 0; t < nPart; t++) { // particle 0 is mother particle
    masses[t] = physConst->FindParticle(finalS.at(t)).GetMass();
  }
  event.SetDecay(W, nPart, masses);
};

RootGenerator::RootGenerator(std::shared_ptr<Kinematics> kin, int seed) {
  gRandom = new TRandom3(0);
  if (seed != -1)
    SetSeed(seed);
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
        "RootGenerator::RootGenerator() | More than one "
        "particle in initial State! Currently we can not deal with that!");

  sqrtS = physConst->FindParticle(initialS.at(0)).GetMass();

  masses = new Double_t[nPart];
  TLorentzVector W(0.0, 0.0, 0.0, sqrtS);    //= beam + target;
  for (unsigned int t = 0; t < nPart; t++) { // particle 0 is mother particle
    masses[t] = physConst->FindParticle(finalS.at(t)).GetMass();
  }
  event.SetDecay(W, nPart, masses);
};

RootGenerator *RootGenerator::Clone() { return (new RootGenerator(*this)); }

void RootGenerator::Generate(Event &evt) {
  evt.Clear();
  const double weight = event.Generate();

  for (unsigned int t = 0; t < nPart; t++) {
    TLorentzVector *p = event.GetDecay(t);
    evt.AddParticle(Particle(p->X(), p->Y(), p->Z(), p->E()));
  }
  evt.SetWeight(weight);

#ifndef _NDEBUG
  ComPWA::FourMomentum pFour;
  for (int i = 0; i < evt.GetNParticles(); i++)
    pFour += evt.GetParticle(i).GetFourMomentum();
  if (pFour.GetInvMass() != sqrtS)
    if (!ComPWA::equal(pFour.GetInvMass(), sqrtS, 100)) {
      LOG(error) << pFour.GetInvMass() << " - " << sqrtS << " = "
                 << pFour.GetInvMass() - sqrtS;
      throw std::runtime_error(
          "RootGenerator::Generate() | Invariant mass of "
          "all generate particles does not sum up to the mass of the decaying "
          "particle. "
          "The difference in larger than two times the numerical precision.");
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
} /* namespace Tools */
} /* namespace ComPWA */
