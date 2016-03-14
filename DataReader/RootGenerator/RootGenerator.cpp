//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------
#include <stdexcept>

#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"

#include "DataReader/RootGenerator/RootGenerator.hpp"

namespace ComPWA {
namespace DataReader {
namespace RootGenerator {

RootGenerator::RootGenerator(int seed) {
  gRandom = new TRandom3(0);
  if (seed != -1)
    setSeed(seed);
  Kinematics* kin = Kinematics::instance();
  nPart = kin->getNumberOfParticles();
  if (nPart < 2)
    throw std::runtime_error(
        "RootGenerator::RootGenerator() | one particle is not enough!");
  if (nPart == 2)
    BOOST_LOG_TRIVIAL(info)<< "RootGenerator::RootGenerator() | only 2 particles in the final"
    " state! There are no degrees of freedom!";
  masses = new Double_t[nPart];
  TLorentzVector W(0.0, 0.0, 0.0, kin->getMotherMass());    //= beam + target;
  for (unsigned int t = 0; t < nPart; t++) {    // particle 0 is mother particle
    masses[t] = kin->getMass(t + 1);
  }
  event.SetDecay(W, nPart, masses);
}
;

RootGenerator* RootGenerator::Clone() {
  return (new RootGenerator(*this));
}

void RootGenerator::generate(Event& evt) {
  bool regenerate_event(false);

  do {
    regenerate_event = false;
    const double weight = event.Generate();
    Event tmp(weight, "");

    for (unsigned int t = 0; t < nPart; t++) {
      TLorentzVector* p = event.GetDecay(t);
      if (p->M2() < 0.0) {
        regenerate_event = true;
        break;
      }
      tmp.addParticle(Particle(p->X(), p->Y(), p->Z(), p->E()));

      if (tmp.getParticle(t).getMassSquare() < 0.0) {
        regenerate_event = true;
        break;
      }
    }
    evt = tmp;
  }
  while (regenerate_event);

  //dataPoint p(tmp); std::cout<<"-"<<p.getVal(0)<<std::endl;
}

void RootGenerator::setSeed(unsigned int seed) {
  gRandom->SetSeed(seed);
}

unsigned int RootGenerator::getSeed() {
  return gRandom->GetSeed();
}
double RootGenerator::getUniform() {
  return gRandom->Uniform(0, 1);
}

void UniformTwoBodyGenerator::generate(Event& evt) {
  double s = RootGenerator::getUniform() * (maxSq - minSq) + minSq;
  TLorentzVector W(0.0, 0.0, 0.0, sqrt(s));    //= beam + target;
//	std::cout<<"generate at s="<<s<<std::endl;
  RootGenerator::getGenerator()->SetDecay(W, nPart, masses);
  RootGenerator::generate(evt);
}

} /* namespace RootGenerator */
} /* namespace DataReader */
} /* namespace ComPWA */
