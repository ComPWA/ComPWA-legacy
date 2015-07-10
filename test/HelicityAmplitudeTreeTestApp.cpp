#include <fstream>

#include "Physics/HelicityAmplitude/HelicityDecayConfiguration.hpp"
#include "Physics/HelicityAmplitude/HelicityDecayXMLConfigReader.hpp"
#include "Physics/HelicityAmplitude/HelicityDecayTreeFactory.hpp"
#include "Physics/HelicityAmplitude/HelicityAmplitudeTreeFactory.hpp"

int main(int argc, char **argv) {
  std::cout << "  ComPWA Copyright (C) 2013  Stefan Pflueger " << std::endl;
  std::cout
      << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt"
      << std::endl;
  std::cout << std::endl;

  HelicityFormalism::HelicityDecayConfiguration decay_configuration;
  HelicityFormalism::HelicityDecayXMLConfigReader xml_reader(
      decay_configuration);
  xml_reader.readConfig("Physics/HelicityAmplitude/JPSI_ypipi.xml");

  HelicityFormalism::HelicityDecayTreeFactory decay_tree_factory(
      decay_configuration);

  std::vector<std::vector<HelicityFormalism::HelicityDecayTree> > decay_trees =
      decay_tree_factory.createDecayTrees();

  HelicityFormalism::HelicityAmplitudeTreeFactory amp_tree_factory;

  std::vector<std::vector<HelicityFormalism::HelicityDecayTree> >::iterator decay_topology;
  for (decay_topology = decay_trees.begin();
      decay_topology != decay_trees.end(); ++decay_topology) {
    std::vector<HelicityFormalism::HelicityDecayTree>::iterator decay_tree;
    for (decay_tree = decay_topology->begin();
        decay_tree != decay_topology->end(); ++decay_tree) {
      if(!decay_tree->hasCycles()) {
        amp_tree_factory.generateAmplitudeTree(*decay_tree);
      }
    }
  }

  return 0;
}

