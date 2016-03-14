#include <fstream>

#include "Core/PhysConst.hpp"

#include "Physics/DecayTree/DecayGenerator.hpp"
#include "Physics/DecayTree/DecayGeneratorFacade.hpp"
#include "Physics/DecayTree/DecayXMLConfigReader.hpp"

int main(int argc, char **argv) {
  std::cout << "  ComPWA Copyright (C) 2013  Stefan Pflueger " << std::endl;
  std::cout
      << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt"
      << std::endl;
  std::cout << std::endl;

  ComPWA::Physics::DecayTree::DecayGenerator decay_generator;
  // initialize

  ComPWA::Physics::DecayTree::DecayGeneratorFacade decay_generator_facade(
      decay_generator);

  decay_generator_facade.setAllowedSpinQuantumNumbers(
      ComPWA::QuantumNumbers::SPIN, { 0, 1, 2}, 1, {
          ComPWA::QuantumNumbers::ORBITAL_ANGULAR_MOMENTUM });
  decay_generator_facade.setAllowedSpinQuantumNumbers(
      ComPWA::QuantumNumbers::ORBITAL_ANGULAR_MOMENTUM, { 0, 1, 2 }, 1, { },
      ComPWA::Physics::DecayTree::QuantumNumberTypes::COMPOSITE_PARTICLE_BASED);
  decay_generator_facade.setAllowedSpinQuantumNumbers(
      ComPWA::QuantumNumbers::ISOSPIN, { 0, 1 }, 1, { });
  decay_generator_facade.setAllowedIntQuantumNumbers(
      ComPWA::QuantumNumbers::CHARGE, { -1, 0, 1 }, { });
  decay_generator_facade.setAllowedIntQuantumNumbers(
      ComPWA::QuantumNumbers::PARITY, { -1, 1 }, {
          ComPWA::QuantumNumbers::ORBITAL_ANGULAR_MOMENTUM });
  decay_generator_facade.setAllowedIntQuantumNumbers(
      ComPWA::QuantumNumbers::CPARITY, { -1, 1 }, {
          ComPWA::QuantumNumbers::CHARGE });

  decay_generator_facade.setConservedQuantumNumbers( {
      ComPWA::QuantumNumbers::CHARGE, ComPWA::QuantumNumbers::PARITY,
      ComPWA::QuantumNumbers::CPARITY, ComPWA::QuantumNumbers::SPIN,
      ComPWA::QuantumNumbers::ORBITAL_ANGULAR_MOMENTUM });

  ComPWA::Physics::DecayTree::IFParticleInfo if_particle =
      decay_generator.createIFParticleInfo("gamma");
  decay_generator.addFinalStateParticles(if_particle);
  if_particle = decay_generator.createIFParticleInfo("pi0");
  decay_generator.addFinalStateParticles(if_particle);
  decay_generator.addFinalStateParticles(if_particle);

  if_particle = decay_generator.createIFParticleInfo("jpsi");
  decay_generator.setTopNodeState(if_particle);

  ComPWA::Physics::DecayTree::DecayConfiguration decay_config =
      decay_generator.createDecayConfiguration();

  ComPWA::Physics::DecayTree::DecayXMLConfigReader xml_config_io(decay_config);
  xml_config_io.writeConfig("test.xml");

  return 0;
}

