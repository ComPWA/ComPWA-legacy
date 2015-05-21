//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#include "HelicityDecayTreeFactory.hpp"

namespace HelicityFormalism {

HelicityDecayTreeFactory::HelicityDecayTreeFactory(
    const HelicityDecayConfiguration& decay_configuration) :
    decay_configuration_(decay_configuration) {
}

HelicityDecayTreeFactory::~HelicityDecayTreeFactory() {
}

bool HelicityDecayTreeFactory::isNodeADaughter(
    DecayTopology::const_iterator& node,
    const std::vector<unsigned int>& list_of_daughters) const {
  for (unsigned int daugther_index = 0;
      daugther_index != list_of_daughters.size(); ++daugther_index) {
    if (list_of_daughters[daugther_index] == node->first)
      return true;
  }
  return false;
}

bool HelicityDecayTreeFactory::isNodeADaughterInTopology(
    DecayTopology::const_iterator& node,
    const DecayTopology& decay_topology) const {
  bool is_daughter(false);
  DecayTopology::const_iterator other_decay_node;
  // check if this index appears in any of the daughter lists if not this is our top node
  for (other_decay_node = decay_topology.begin();
      other_decay_node != decay_topology.end(); ++other_decay_node) {
    is_daughter = isNodeADaughter(node, other_decay_node->second.first);
    if (is_daughter)
      break;
    is_daughter = isNodeADaughter(node, other_decay_node->second.second);
    if (is_daughter)
      break;
  }
  return is_daughter;
}

DecayTopology::const_iterator HelicityDecayTreeFactory::determineTopNode(
    const DecayTopology& decay_topology) const {
  std::vector<DecayTopology::const_iterator> candidates;
  // loop through the topology
  DecayTopology::const_iterator decay_node;

  for (decay_node = decay_topology.begin(); decay_node != decay_topology.end();
      ++decay_node) {
    // if this node is not a daughter in any other decay
    // then add this to the list of candidates
    if (!isNodeADaughterInTopology(decay_node, decay_topology)) {
      candidates.push_back(decay_node);
    }
  }
  // if we have more than 1 candidate throw an error...
  if (candidates.size() > 1)
    throw std::runtime_error(
        "The decay topology has more than one top node, meaning the topology is corrupted. Please fix the configuration file and rerun!");
  if (candidates.size() == 0)
    throw std::runtime_error(
        "The decay topology does not have any top node, meaning the topology is corrupted. Please fix the configuration file and rerun!");
  // otherwise return here
  return candidates[0];
}

bool HelicityDecayTreeFactory::canDecayTreesGrow(
    const std::vector<HelicityDecayTree>& decay_trees,
    const DecayTopology& decay_topology) const {
  bool has_nodes_to_grow(false);

  // loop over the decay trees
  std::vector<HelicityDecayTree>::const_iterator decay_tree;
  for (decay_tree = decay_trees.begin(); decay_tree != decay_trees.end();
      ++decay_tree) {
    // check if there are decays available for the current leaves

    std::vector<ParticleState> list_of_leaves = decay_tree->getLowestLeaves();

    for (unsigned int leaf_index = 0; leaf_index < list_of_leaves.size();
        ++leaf_index) {
      std::vector<ParticleState>::const_iterator particle = std::find(
          decay_configuration_.particles_.begin(),
          decay_configuration_.particles_.end(), list_of_leaves[leaf_index]);
      DecayTopology::const_iterator found_decay = decay_topology.find(
          particle->particle_id_);
      if (found_decay != decay_topology.end()) {
        has_nodes_to_grow = true;
      }
    }
  }
  return has_nodes_to_grow;
}

std::vector<HelicityDecayTree> HelicityDecayTreeFactory::growSingleLeafOnDecayTrees(
    const std::vector<HelicityDecayTree>& decay_trees,
    DecayTopology::const_iterator& single_decay) const {
  std::vector<HelicityDecayTree> new_decay_trees;

  // loop over the decay trees
  std::vector<HelicityDecayTree>::const_iterator decay_tree;
  for (decay_tree = decay_trees.begin(); decay_tree != decay_trees.end();
      ++decay_tree) {
    const ParticleState& particle =
        decay_configuration_.particles_[single_decay->first];
    for (unsigned int daugther_1_index = 0;
        daugther_1_index != single_decay->second.first.size();
        ++daugther_1_index) {
      for (unsigned int daugther_2_index = 0;
          daugther_2_index != single_decay->second.second.size();
          ++daugther_2_index) {
        HelicityDecayTree temp_tree(*decay_tree);
        temp_tree.createDecay(
            decay_configuration_.particles_[single_decay->first],
            std::make_pair(
                decay_configuration_.particles_[single_decay->second.first[daugther_1_index]],
                decay_configuration_.particles_[single_decay->second.second[daugther_2_index]]));
        new_decay_trees.push_back(temp_tree);
      }
    }
  }
  return new_decay_trees;
}

std::vector<HelicityDecayTree> HelicityDecayTreeFactory::growNextDecayTreeLayer(
    const std::vector<HelicityDecayTree>& decay_trees,
    const DecayTopology& decay_topology) const {
  std::vector<HelicityDecayTree> new_decay_trees;

// loop over the decay trees
  std::vector<HelicityDecayTree>::const_iterator decay_tree;
  for (decay_tree = decay_trees.begin(); decay_tree != decay_trees.end();
      ++decay_tree) {
    // check if there are decays available for the current leaves
    std::vector<ParticleState> list_of_leaves = decay_tree->getLowestLeaves();
    // TODO: this is bad practice but works right now... clear the list
    std::vector<HelicityDecayTree> temp_trees;
    temp_trees.push_back(*decay_tree);
    temp_trees[0].clearCurrentGrownNodes();

    for (unsigned int leaf_index = 0; leaf_index < list_of_leaves.size();
        ++leaf_index) {
      std::vector<ParticleState>::const_iterator particle = std::find(
          decay_configuration_.particles_.begin(),
          decay_configuration_.particles_.end(), list_of_leaves[leaf_index]);
      unsigned int particle_index = std::distance(
          decay_configuration_.particles_.begin(), particle);
      DecayTopology::const_iterator found_decay = decay_topology.find(
          particle_index);
      if (found_decay != decay_topology.end()) {
        temp_trees = growSingleLeafOnDecayTrees(temp_trees, found_decay);
      }
    }
    new_decay_trees.insert(new_decay_trees.end(), temp_trees.begin(),
        temp_trees.end());
  }
  return new_decay_trees;
}

std::vector<HelicityDecayTree> HelicityDecayTreeFactory::createDecayTreeSeedList(
    DecayTopology::const_iterator & top_node_iter) const {
  std::vector<HelicityDecayTree> decay_trees;
  for (unsigned int daugther_1_index = 0;
      daugther_1_index != top_node_iter->second.first.size();
      ++daugther_1_index) {
    for (unsigned int daugther_2_index = 0;
        daugther_2_index != top_node_iter->second.second.size();
        ++daugther_2_index) {
      HelicityDecayTree temp_tree;
      temp_tree.createDecay(
          decay_configuration_.particles_[top_node_iter->first],
          std::make_pair(
              decay_configuration_.particles_[top_node_iter->second.first[daugther_1_index]],
              decay_configuration_.particles_[top_node_iter->second.second[daugther_2_index]]));
      decay_trees.push_back(temp_tree);
    }
  }
  return decay_trees;
}

std::vector<std::vector<HelicityDecayTree> > HelicityDecayTreeFactory::createDecayTrees() const {
  std::vector<std::vector<HelicityDecayTree> > decay_trees;

// loop over the topologies
  std::vector<DecayTopology>::const_iterator decay_topology_iter;
  for (decay_topology_iter = decay_configuration_.decay_topologies_.begin();
      decay_topology_iter != decay_configuration_.decay_topologies_.end();
      ++decay_topology_iter) {

    // determine top node in the current topology
    DecayTopology::const_iterator top_node = determineTopNode(
        *decay_topology_iter);

    // then start the recursive tree builder algorithm with that node
    // set seeds
    std::vector<HelicityDecayTree> temp_decay_tree_list =
        createDecayTreeSeedList(top_node);

    // grow trees as much as possible
    while (canDecayTreesGrow(temp_decay_tree_list, *decay_topology_iter)) {
      temp_decay_tree_list = growNextDecayTreeLayer(temp_decay_tree_list,
          *decay_topology_iter);;
    }

    // and simply add the list of constructed trees to the full decay tree list
    decay_trees.push_back(temp_decay_tree_list);
  }

  return decay_trees;
}

} /* namespace HelicityFormalism */
