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

#include "Core/Utility.hpp"

#include "Physics/DecayTree/DecayTreeFactory.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayTree {

DecayTreeFactory::DecayTreeFactory(
    const DecayConfiguration& decay_configuration) :
    decay_configuration_(decay_configuration) {
}

DecayTreeFactory::~DecayTreeFactory() {
}

std::vector<DecayTree> DecayTreeFactory::createDecayTrees() const {
  std::vector<DecayTree> decay_trees;

  // loop over the decay trees of the configuration
  std::vector<ParticleIndexDecayTree>::const_iterator decay_tree_iter;
  for (decay_tree_iter = decay_configuration_.concrete_decay_trees_.begin();
      decay_tree_iter != decay_configuration_.concrete_decay_trees_.end();
      ++decay_tree_iter) {

    // determine top node in the current topology
    ParticleIndexDecayTree::const_iterator top_node =
        decay_configuration_.determineTopNode(*decay_tree_iter);

    // then start the recursive tree builder algorithm with that node
    // set seeds
    DecayTree decay_tree = createSeedDecayTree(top_node);

    // grow tree as much as possible
    while (canDecayTreeGrow(decay_tree, *decay_tree_iter)) {
      growNextDecayTreeLayer(decay_tree, *decay_tree_iter);
    }

    if (!isDecayTreeFaulty(decay_tree)) {
      // and simply add the list of constructed trees to the full decay tree list
      decay_trees.push_back(decay_tree);
    }
  }

  removeUndistinguishableCombinations(decay_trees);

  return decay_trees;
}

DecayTree DecayTreeFactory::createSeedDecayTree(
    ParticleIndexDecayTree::const_iterator & top_node_iter) const {
  DecayTree decay_tree;

  growSingleLeafOnDecayTree(decay_tree, top_node_iter);

  return decay_tree;
}

void DecayTreeFactory::growSingleLeafOnDecayTree(DecayTree& decay_tree,
    ParticleIndexDecayTree::const_iterator& single_decay) const {

  std::vector<ParticleStateInfo> decay_daughters(
      generateParticleStateInfoList(single_decay->second.particle_indices_));

  DecayNode mother_node;
  mother_node.state_info_ =
      decay_configuration_.particles_[single_decay->first];
  mother_node.strength_and_phase_ =
      single_decay->second.decay_strength_info_and_phase_;
  decay_tree.createDecay(mother_node, decay_daughters);
}

std::vector<ParticleStateInfo> DecayTreeFactory::generateParticleStateInfoList(
    const std::vector<unsigned int>& particle_index_list) const {
  std::vector<ParticleStateInfo> particle_list;

  for (unsigned int i = 0; i < particle_index_list.size(); ++i) {
    particle_list.push_back(
        decay_configuration_.particles_[particle_index_list[i]]);
  }
  return particle_list;
}

bool DecayTreeFactory::canDecayTreeGrow(const DecayTree& decay_tree,
    const ParticleIndexDecayTree& decay_topology) const {
  bool has_nodes_to_grow(false);

  // check if there are decays available for the current leaves
  std::vector<DecayNode> list_of_leaves = decay_tree.getLowestLeaves();

  for (unsigned int leaf_index = 0; leaf_index < list_of_leaves.size();
      ++leaf_index) {
    std::vector<ParticleStateInfo>::const_iterator particle = std::find(
        decay_configuration_.particles_.begin(),
        decay_configuration_.particles_.end(),
        list_of_leaves[leaf_index].state_info_);
    unsigned int particle_index = particle
        - decay_configuration_.particles_.begin();
    ParticleIndexDecayTree::const_iterator found_decay = decay_topology.find(
        particle_index);
    if (found_decay != decay_topology.end()) {
      has_nodes_to_grow = true;
    }
  }
  return has_nodes_to_grow;
}

void DecayTreeFactory::growNextDecayTreeLayer(DecayTree& decay_tree,
    const ParticleIndexDecayTree& decay_topology) const {

  // check if there are decays available for the current leaves
  std::vector<DecayNode> list_of_leaves = decay_tree.getLowestLeaves();
  decay_tree.clearCurrentGrownNodes();

  for (unsigned int leaf_index = 0; leaf_index < list_of_leaves.size();
      ++leaf_index) {
    std::vector<ParticleStateInfo>::const_iterator particle = std::find(
        decay_configuration_.particles_.begin(),
        decay_configuration_.particles_.end(),
        list_of_leaves[leaf_index].state_info_);
    unsigned int particle_index = particle
        - decay_configuration_.particles_.begin();

    ParticleIndexDecayTree::const_iterator found_decay = decay_topology.find(
        particle_index);
    if (found_decay != decay_topology.end()) {
      growSingleLeafOnDecayTree(decay_tree, found_decay);
    }
  }
}

bool DecayTreeFactory::isDecayTreeFaulty(const DecayTree& decay_tree) const {
  if (decay_tree.hasCycles() || decay_tree.isDisconnected()) {
    return true;
  }
  return false;
}

void DecayTreeFactory::removeUndistinguishableCombinations(
    std::vector<DecayTree>& decay_tree_list) const {
  std::map<std::vector<std::string>, DecayTree> unique_combinations;
  for (unsigned int i = 0; i < decay_tree_list.size(); ++i) {
    std::vector<std::string> particle_id_combinations = convertToUniqueNameList(
        decay_tree_list[i]);

    if (unique_combinations.find(particle_id_combinations)
        == unique_combinations.end()) {
      unique_combinations[particle_id_combinations] = decay_tree_list[i];
    }
  }

  decay_tree_list.clear();
  for (auto unique_combination_iter = unique_combinations.begin();
      unique_combination_iter != unique_combinations.end();
      ++unique_combination_iter) {
    decay_tree_list.push_back(unique_combination_iter->second);
  }
}

std::vector<std::string> DecayTreeFactory::convertToUniqueNameList(
    const DecayTree& decay_tree) const {
  std::vector<std::string> name_list;

  std::pair<boost::graph_traits<HelicityTree>::vertex_iterator,
      boost::graph_traits<HelicityTree>::vertex_iterator> vp;
  for (vp = vertices(decay_tree.getHelicityDecayTree()); vp.first != vp.second;
      ++vp.first) {
    name_list.push_back(
        decay_tree.getHelicityDecayTree()[*vp.first].state_info_.pid_information_.name_);
  }
  std::sort(name_list.begin(), name_list.end());
  return name_list;
}

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */
