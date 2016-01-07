//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
//----------------------------------------------------------------------------------

#include <fstream>

#include "Physics/DecayTree/DecayConfiguration.hpp"
#include "Physics/DecayTree/DecayXMLConfigReader.hpp"
#include "Physics/DecayTree/DecayTreeFactory.hpp"

int main(int argc, char **argv) {
  std::cout << "  ComPWA Copyright (C) 2013  Stefan Pflueger " << std::endl;
  std::cout
      << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt"
      << std::endl;
  std::cout << std::endl;

  std::string input_config_file("Physics/HelicityAmplitude/JPSI_ypipi.xml");
  std::string output_file("graph.dot");

  ComPWA::DecayTree::DecayConfiguration decay_configuration;
  ComPWA::DecayTree::DecayXMLConfigReader xml_reader(decay_configuration);
  xml_reader.readConfig(input_config_file);

  ComPWA::DecayTree::DecayTreeFactory decay_tree_factory(decay_configuration);

  std::vector<ComPWA::DecayTree::DecayTree> decay_trees =
      decay_tree_factory.createDecayTrees();

  std::cout << "created " << decay_trees.size() << " decay trees from "
      << input_config_file << " config file!" << std::endl;
  std::cout << "printing to " << output_file << " output file" << std::endl;

  std::ofstream dot(output_file);

  std::vector<ComPWA::DecayTree::DecayTree>::iterator decay_tree;
  for (decay_tree = decay_trees.begin(); decay_tree != decay_trees.end();
      ++decay_tree) {
    decay_tree->print(dot);
  }

  return 0;
}

