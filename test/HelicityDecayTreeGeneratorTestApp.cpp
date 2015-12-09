#include <fstream>

#include "Physics/HelicityAmplitude/DecayGenerator.hpp"

int main(int argc, char **argv) {
  std::cout << "  ComPWA Copyright (C) 2013  Stefan Pflueger " << std::endl;
  std::cout
      << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt"
      << std::endl;
  std::cout << std::endl;

  HelicityFormalism::DecayGenerator decay_generator;
  // initialize
  decay_generator.addFinalStateParticles("gamma");
//  decay_generator.addFinalStateParticles("gamma");
//  decay_generator.addFinalStateParticles("pi0");
  decay_generator.addFinalStateParticles("pi0");
  decay_generator.addFinalStateParticles("pi0");

  decay_generator.addIntermediateStateParticles("f0_980");
  decay_generator.addIntermediateStateParticles("omega");

  decay_generator.setTopNodeState("jpsi");

  decay_generator.generate();

  return 0;
}

