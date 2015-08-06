#include <fstream>

#include "Physics/HelicityAmplitude/DecayConfiguration.hpp"
#include "Physics/HelicityAmplitude/DecayXMLConfigReader.hpp"
#include "Physics/HelicityAmplitude/DecayTreeFactory.hpp"
#include "Physics/HelicityAmplitude/TopologyAmplitudeFactory.hpp"

int main(int argc, char **argv) {
  std::cout << "  ComPWA Copyright (C) 2013  Stefan Pflueger " << std::endl;
  std::cout
      << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt"
      << std::endl;
  std::cout << std::endl;

  HelicityFormalism::DecayConfiguration decay_configuration;
  HelicityFormalism::DecayXMLConfigReader xml_reader(
      decay_configuration);
  xml_reader.readConfig("Physics/HelicityAmplitude/JPSI_ypipi.xml");

  HelicityFormalism::DecayTreeFactory decay_tree_factory(
      decay_configuration);

  std::vector<HelicityFormalism::DecayTree> decay_trees =
      decay_tree_factory.createDecayTrees();

  HelicityFormalism::TopologyAmplitudeFactory amp_tree_factory;

  amp_tree_factory.generateTopologyAmplitudes(decay_trees);

  return 0;
}

