#include <fstream>

#include "Core/PhysConst.hpp"

#include "Physics/DecayTree/DecayGenerator.hpp"
#include "Physics/DecayTree/DecayGeneratorFacade.hpp"

int main(int argc, char **argv) {
  std::cout << "  ComPWA Copyright (C) 2013  Stefan Pflueger " << std::endl;
  std::cout
      << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt"
      << std::endl;
  std::cout << std::endl;

  ComPWA::DecayTree::DecayGenerator decay_generator;
  // initialize

  ComPWA::DecayTree::DecayGeneratorFacade decay_generator_facade(
      decay_generator);

  decay_generator_facade.setAllowedSpinQuantumNumbers(
      ComPWA::QuantumNumbers::SPIN, { 0, 1, 2 }, 1);
  decay_generator_facade.setAllowedSpinQuantumNumbers(
      ComPWA::QuantumNumbers::ISOSPIN, { 0, 1 }, 1);
  decay_generator_facade.setAllowedIntQuantumNumbers(
      ComPWA::QuantumNumbers::CHARGE, { -1, 0, 1 });
  decay_generator_facade.setAllowedIntQuantumNumbers(
      ComPWA::QuantumNumbers::PARITY, { -1, 1 });
  decay_generator_facade.setAllowedIntQuantumNumbers(
      ComPWA::QuantumNumbers::CPARITY, { -1, 1 });

  ComPWA::DecayTree::IFParticleInfo if_particle =
      decay_generator.createIFParticleInfo("gamma");
  decay_generator.addFinalStateParticles(if_particle);
//  decay_generator.addFinalStateParticles("gamma");
//  decay_generator.addFinalStateParticles("pi0");
  if_particle = decay_generator.createIFParticleInfo("pi0");
  decay_generator.addFinalStateParticles(if_particle);
  decay_generator.addFinalStateParticles(if_particle);

  if_particle = decay_generator.createIFParticleInfo("jpsi");
  decay_generator.setTopNodeState(if_particle);

  ComPWA::DecayTree::DecayConfiguration decay_config = decay_generator.createDecayConfiguration();

  return 0;
}

