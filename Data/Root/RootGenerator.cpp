// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Data/Root/RootGenerator.hpp"

#include "Core/Event.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "Core/Random.hpp"
#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"

namespace ComPWA {
namespace Data {
namespace Root {

RootUniformRealGenerator::RootUniformRealGenerator(int seed)
    : RandomGenerator(seed), Seed(seed) {}

double RootUniformRealGenerator::operator()() { return RandomGenerator.Rndm(); }

int RootUniformRealGenerator::getSeed() const { return Seed; }

void RootUniformRealGenerator::setSeed(int seed) {
  Seed = seed;
  RandomGenerator.SetSeed(seed);
}

RootGenerator::RootGenerator(const ComPWA::FourMomentum &CMSP4_,
                             const std::vector<double> &FinalStateMasses_,
                             const std::vector<ComPWA::pid> &FinalStatePIDs_)
    : CMSP4(CMSP4_), FinalStateMasses(FinalStateMasses_),
      FinalStatePIDs(FinalStatePIDs_), CMSBoostVector(0.0, 0.0, 0.0) {
  unsigned int nPart = FinalStateMasses.size();
  if (nPart < 2)
    throw std::runtime_error(
        "RootGenerator::RootGenerator() | one particle is not enough!");
  if (nPart == 2)
    LOG(INFO)
        << "RootGenerator::RootGenerator() | only 2 particles in the final"
           " state! There are no degrees of freedom!";

  init();
}

RootGenerator::RootGenerator(const ComPWA::ParticleList &PartL,
                             std::vector<pid> InitialS, std::vector<pid> FinalS)
    : RootGenerator(
          [&]() {
            if (InitialS.size() != 1)
              throw std::runtime_error(
                  "RootGenerator::RootGenerator() | More than one "
                  "particle in initial State! You need to specify a cms "
                  "four-momentum");
            return ComPWA::FourMomentum(
                0.0, 0.0, 0.0,
                findParticle(PartL, InitialS.at(0)).getMass().Value);
          }(),
          [&]() {
            std::vector<double> fsm;
            for (auto ParticlePid : FinalS) { // particle 0 is mother particle
              fsm.push_back(findParticle(PartL, ParticlePid).getMass().Value);
            }
            return fsm;
          }(),
          FinalS) {}

RootGenerator::RootGenerator(
    const Physics::ParticleStateTransitionKinematicsInfo &KinematicsInfo)
    : RootGenerator(KinematicsInfo.getInitialStateFourMomentum(),
                    KinematicsInfo.getFinalStateMasses(),
                    KinematicsInfo.getFinalStatePIDs()) {}

void RootGenerator::init() {
  CMSEnergyMinusMasses = CMSP4.invariantMass();
  for (double fsmass : FinalStateMasses) {
    CMSEnergyMinusMasses -= fsmass;
  }

  if (CMSEnergyMinusMasses <= 0)
    throw std::runtime_error(
        "RootGenerator::init(): not enough energy for this decay");

  //
  //------> the max weight depends on opt:
  //   opt == "Fermi"  --> fermi energy dependence for cross section
  //   else            --> constant cross section as function of TECM (default)
  //
  /*if (strcasecmp(opt, "fermi") == 0) {
    // ffq[] = pi * (2*pi)**(FNt-2) / (FNt-2)!
    Double_t ffq[] = {0,        3.141592, 19.73921, 62.01255, 129.8788,
                      204.0131, 256.3704, 268.4705, 240.9780, 189.2637,
                      132.1308, 83.0202,  47.4210,  24.8295,  12.0006,
                      5.3858,   2.2560,   0.8859};
    fWtMax = TMath::Power(fTeCmTm, fNt - 2) * ffq[fNt - 1] / P.Mag();

  } else {*/
  double emmax = CMSEnergyMinusMasses + FinalStateMasses[0];
  double emmin = 0.0;
  double wtmax = 1.0;
  for (unsigned int n = 1; n < FinalStateMasses.size(); ++n) {
    emmin += FinalStateMasses[n - 1];
    emmax += FinalStateMasses[n];
    wtmax *= PDK(emmax, emmin, FinalStateMasses[n]);
  }
  MaximumWeight = 1.0 / wtmax;

  TLorentzVector W(CMSP4.px(), CMSP4.py(), CMSP4.pz(), CMSP4.e());
  if (W.Beta()) {
    double w = W.Beta() / W.Rho();
    CMSBoostVector.SetXYZ(W(0) * w, W(1) * w, W(2) * w);
  }
}

ComPWA::EventCollection
RootGenerator::generate(unsigned int NumberOfEvents,
                        UniformRealNumberGenerator &RandomGenerator) const {
  EventCollection GeneratedPhsp{FinalStatePIDs};

  for (unsigned int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {

    size_t NumberOfFinalStateParticles(FinalStateMasses.size());
    std::vector<double> OrderedRandomNumbers;
    OrderedRandomNumbers.reserve(NumberOfFinalStateParticles);
    OrderedRandomNumbers.push_back(0.0);

    if (NumberOfFinalStateParticles > 2) {
      for (unsigned int n = 1; n < NumberOfFinalStateParticles - 1; ++n)
        OrderedRandomNumbers.push_back(RandomGenerator()); // N-2 random numbers
      std::sort(OrderedRandomNumbers.begin(), OrderedRandomNumbers.end());
    }
    OrderedRandomNumbers.push_back(1.0);

    double invMas[NumberOfFinalStateParticles], sum = 0.0;
    for (unsigned int n = 0; n < NumberOfFinalStateParticles; ++n) {
      sum += FinalStateMasses[n];
      invMas[n] = OrderedRandomNumbers[n] * CMSEnergyMinusMasses + sum;
    }

    // compute the weight of the current event
    double weight = MaximumWeight;
    double pd[NumberOfFinalStateParticles];
    for (unsigned int n = 0; n < NumberOfFinalStateParticles - 1; ++n) {
      pd[n] = PDK(invMas[n + 1], invMas[n], FinalStateMasses[n + 1]);
      weight *= pd[n];
    }

    std::vector<TLorentzVector> FinalStateLorentzVectors(
        FinalStateMasses.size());

    // complete specification of event (Raubold-Lynch method)
    FinalStateLorentzVectors[0].SetPxPyPzE(
        0.0, pd[0], 0.0,
        std::sqrt(std::pow(pd[0], 2) + std::pow(FinalStateMasses[0], 2)));

    unsigned int i(1);
    while (true) {
      FinalStateLorentzVectors[i].SetPxPyPzE(
          0.0, -pd[i - 1], 0.0,
          std::sqrt(std::pow(pd[i - 1], 2) + std::pow(FinalStateMasses[i], 2)));

      double cZ = 2.0 * RandomGenerator() - 1.0;
      double sZ = std::sqrt(1.0 - std::pow(cZ, 2));
      double angY = 2.0 * TMath::Pi() * RandomGenerator();
      double cY = std::cos(angY);
      double sY = std::sin(angY);
      for (unsigned int j = 0; j <= i; ++j) {
        TLorentzVector &v(FinalStateLorentzVectors[j]);
        double x = v.Px();
        double y = v.Py();
        double z = v.Pz();
        // rotation around Z and Y
        v.SetPx(cY * (cZ * x - sZ * y) - sY * z);
        v.SetPy(sZ * x + cZ * y);
        v.SetPz(sY * (cZ * x - sZ * y) + cY * z);
      }
      if (i == NumberOfFinalStateParticles - 1)
        break;
      // double beta = 1.0 / std::sqrt(1.0 + std::pow(invMas[i] / pd[i], 2));
      double beta2 = 1.0 / (1.0 + std::pow(invMas[i] / pd[i], 2));
      for (unsigned int j = 0; j <= i; ++j) {
        // FinalStateLorentzVectors[j].Boost(0.0, beta, 0.0);
        BoostAlongY(FinalStateLorentzVectors[j], beta2);
      }
      ++i;
    }

    std::vector<FourMomentum> FourMomenta;
    // final boost of all particles
    for (unsigned int n = 0; n < NumberOfFinalStateParticles; ++n) {
      FinalStateLorentzVectors[n].Boost(CMSBoostVector);
      FourMomenta.push_back(FourMomentum(
          FinalStateLorentzVectors[n].X(), FinalStateLorentzVectors[n].Y(),
          FinalStateLorentzVectors[n].Z(), FinalStateLorentzVectors[n].E()));
    }
    GeneratedPhsp.Events.push_back(ComPWA::Event{FourMomenta, 1.0});
  }
  return GeneratedPhsp;
}

void RootGenerator::BoostAlongY(TLorentzVector &vec,
                                double beta_squared) const {
  // Boost this Lorentz vector
  double y(vec.Y());
  double t(vec.T());

  double gamma = 1.0 / std::sqrt(1.0 - beta_squared);
  double betagamma = 1.0 / std::sqrt((1.0 - beta_squared) / beta_squared);

  vec.SetY(gamma * y + betagamma * t);
  vec.SetT(gamma * t + betagamma * y);
}

double RootGenerator::PDK(double a, double b, double c) const {
  double x = (a - b - c) * (a + b + c) * (a - b + c) * (a + b - c);
  return std::sqrt(x) / (2.0 * a);
}

} // namespace Root
} // namespace Data
} // namespace ComPWA
