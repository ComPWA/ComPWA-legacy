#include <fstream>

#include "Physics/DecayTree/DecayGenerator.hpp"

int main(int argc, char **argv) {
  std::cout << "  ComPWA Copyright (C) 2013  Stefan Pflueger " << std::endl;
  std::cout
      << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt"
      << std::endl;
  std::cout << std::endl;

  ComPWA::DecayTree::DecayGenerator decay_generator;
  // initialize
  ComPWA::DecayTree::IFParticleInfo if_particle = decay_generator.createIFParticleInfo("gamma");
  decay_generator.addFinalStateParticles(if_particle);
//  decay_generator.addFinalStateParticles("gamma");
//  decay_generator.addFinalStateParticles("pi0");
  if_particle = decay_generator.createIFParticleInfo("pi0");
  decay_generator.addFinalStateParticles(if_particle);
  decay_generator.addFinalStateParticles(if_particle);

  if_particle = decay_generator.createIFParticleInfo("jpsi");
  decay_generator.setTopNodeState(if_particle);

  decay_generator.generate();

  return 0;
}

