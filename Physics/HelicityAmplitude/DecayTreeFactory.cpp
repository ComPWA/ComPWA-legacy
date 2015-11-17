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

#include "DecayTreeFactory.hpp"
#include "ParticleStateDefinitions.hpp"

namespace HelicityFormalism {

DecayTreeFactory::DecayTreeFactory(
    const DecayConfiguration& decay_configuration) :
    decay_configuration_(decay_configuration) {
}

DecayTreeFactory::~DecayTreeFactory() {
}

std::vector<DecayTree> DecayTreeFactory::createDecayTrees() const {
  std::vector<DecayTree> decay_trees;

  // loop over the decay trees of the configuration
  std::vector<ParticleIndexDecayTree>::const_iterator decay_topology_iter;
  for (decay_topology_iter = decay_configuration_.concrete_decay_trees_.begin();
      decay_topology_iter != decay_configuration_.concrete_decay_trees_.end();
      ++decay_topology_iter) {

    // determine top node in the current topology
    ParticleIndexDecayTree::const_iterator top_node = determineTopNode(
        *decay_topology_iter);

    // then start the recursive tree builder algorithm with that node
    // set seeds
    std::vector<DecayTree> temp_decay_tree_list = createDecayTreeSeedList(
        top_node);

    // grow trees as much as possible
    while (canDecayTreesGrow(temp_decay_tree_list, *decay_topology_iter)) {
      temp_decay_tree_list = growNextDecayTreeLayer(temp_decay_tree_list,
          *decay_topology_iter);
    }

    removeFaultyDecayTrees(temp_decay_tree_list);

    // and simply add the list of constructed trees to the full decay tree list
    decay_trees.insert(decay_trees.begin(), temp_decay_tree_list.begin(),
        temp_decay_tree_list.end());
  }

  return decay_trees;
}

ParticleIndexDecayTree::const_iterator DecayTreeFactory::determineTopNode(
    const ParticleIndexDecayTree& decay_topology) const {
  std::vector<ParticleIndexDecayTree::const_iterator> candidates;
  // loop through the topology
  ParticleIndexDecayTree::const_iterator decay_node;

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
  else if (candidates.size() == 0)
    throw std::runtime_error(
        "The decay topology does not have any top node, meaning the topology is corrupted. Please fix the configuration file and rerun!");
  // otherwise return here
  return candidates[0];
}

bool DecayTreeFactory::isNodeADaughterInTopology(
    ParticleIndexDecayTree::const_iterator& node,
    const ParticleIndexDecayTree& decay_topology) const {
  bool is_daughter(false);
  ParticleIndexDecayTree::const_iterator other_decay_node;
  // check if this index appears in any of the daughter lists if not this is our top node
  for (other_decay_node = decay_topology.begin();
      other_decay_node != decay_topology.end(); ++other_decay_node) {
    std::vector<DecayProductsInfo>::const_iterator decay_products_iter;
    for (decay_products_iter = other_decay_node->second.begin();
        decay_products_iter != other_decay_node->second.end();
        ++decay_products_iter) {
      is_daughter = isNodeADaughter(node,
          decay_products_iter->particle_indices_);
      if (is_daughter)
        break;
    }
  }
  return is_daughter;
}

bool DecayTreeFactory::isNodeADaughter(
    ParticleIndexDecayTree::const_iterator& node,
    const std::vector<unsigned int>& list_of_daughters) const {
  for (unsigned int daugther_index = 0;
      daugther_index != list_of_daughters.size(); ++daugther_index) {
    if (list_of_daughters[daugther_index] == node->first)
      return true;
  }
  return false;
}

std::vector<DecayTree> DecayTreeFactory::createDecayTreeSeedList(
    ParticleIndexDecayTree::const_iterator & top_node_iter) const {
  std::vector<DecayTree> decay_trees;

  for (auto decay_products_info = top_node_iter->second.begin();
      decay_products_info != top_node_iter->second.end();
      ++decay_products_info) {
    std::vector<unsigned int> current_particle_index_list();

    std::vector<ParticleStateInfo> decay_daughters(
        generateParticleStateInfoList(decay_products_info->particle_indices_));

    DecayTree temp_tree;
    DecayNode mother_node;
    mother_node.state_info_ =
        decay_configuration_.particles_[top_node_iter->first];
    mother_node.strength_and_phase_ =
        decay_products_info->decay_strength_info_and_phase_;
    temp_tree.createDecay(mother_node, decay_daughters);
    decay_trees.push_back(temp_tree);
  }
  return decay_trees;
}

std::vector<ParticleStateInfo> DecayTreeFactory::generateParticleStateInfoList(
    const std::vector<unsigned int>& particle_index_list) const {
  std::vector<ParticleStateInfo> particle_list;

  for (unsigned int i = 0; i < particle_index_list.size(); ++i) {
    particle_list.push_back(
        decay_configuration_.particles_[particle_index_list[i]]);

    // extendDaughterCombinations(particle_index_lists, initial_particle_list,
    //     decay_daughters_combination_list);
  }
  return particle_list;
}

/*void DecayTreeFactory::extendDaughterCombinations(
 std::vector<std::vector<unsigned int> > remaining_particle_index_lists,
 const std::vector<ParticleStateInfo>& decay_daughters_combination,
 std::vector<std::vector<ParticleStateInfo> >& decay_daughters_combination_list) const {

 if (remaining_particle_index_lists.size() > 0) {
 std::cout<<"wth"<<std::endl;
 std::vector<unsigned int> current_particle_index_list(
 remaining_particle_index_lists.back());
 remaining_particle_index_lists.pop_back();

 for (unsigned int i = 0; i < current_particle_index_list.size(); ++i) {
 std::vector<ParticleStateInfo> current_particle_list(
 decay_daughters_combination);

 const ParticleStateInfo &particle =
 decay_configuration_.particles_[current_particle_index_list[i]];

 if (isParticleValidForCombination(particle, current_particle_list)) {
 current_particle_list.push_back(particle);

 extendDaughterCombinations(remaining_particle_index_lists,
 current_particle_list, decay_daughters_combination_list);
 }
 }
 }
 else {
 decay_daughters_combination_list.push_back(decay_daughters_combination);
 }
 }

 bool DecayTreeFactory::isParticleValidForCombination(
 const ParticleStateInfo& particle,
 const std::vector<ParticleStateInfo>& combination) const {
 // if we find a particle with the same particle id in the list as the
 // reference particle, this will be not valid
 if (std::find_if(combination.begin(), combination.end(),
 ParticleStateIDComparison(particle.id_)) != combination.end()) {
 return false;
 }
 return true;
 }*/

bool DecayTreeFactory::canDecayTreesGrow(
    const std::vector<DecayTree>& decay_trees,
    const ParticleIndexDecayTree& decay_topology) const {
  bool has_nodes_to_grow(false);

// loop over the decay trees
  std::vector<DecayTree>::const_iterator decay_tree;
  for (decay_tree = decay_trees.begin(); decay_tree != decay_trees.end();
      ++decay_tree) {
    // check if there are decays available for the current leaves

    std::vector<DecayNode> list_of_leaves = decay_tree->getLowestLeaves();

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
  }
  return has_nodes_to_grow;
}

std::vector<DecayTree> DecayTreeFactory::growNextDecayTreeLayer(
    const std::vector<DecayTree>& decay_trees,
    const ParticleIndexDecayTree& decay_topology) const {
  std::vector<DecayTree> new_decay_trees;

// loop over the decay trees
  std::vector<DecayTree>::const_iterator decay_tree;
  for (decay_tree = decay_trees.begin(); decay_tree != decay_trees.end();
      ++decay_tree) {
    // check if there are decays available for the current leaves
    std::vector<DecayNode> list_of_leaves = decay_tree->getLowestLeaves();
    // TODO: this is bad practice but works right now... clear the list
    std::vector<DecayTree> temp_trees;
    temp_trees.push_back(*decay_tree);
    temp_trees[0].clearCurrentGrownNodes();

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
        temp_trees = growSingleLeafOnDecayTrees(temp_trees, found_decay);
      }
    }
    new_decay_trees.insert(new_decay_trees.end(), temp_trees.begin(),
        temp_trees.end());
  }
  return new_decay_trees;
}

std::vector<DecayTree> DecayTreeFactory::growSingleLeafOnDecayTrees(
    const std::vector<DecayTree>& decay_trees,
    ParticleIndexDecayTree::const_iterator& single_decay) const {
  std::vector<DecayTree> new_decay_trees;

// loop over the decay trees
  std::vector<DecayTree>::const_iterator decay_tree;
  for (decay_tree = decay_trees.begin(); decay_tree != decay_trees.end();
      ++decay_tree) {
    const ParticleStateInfo& mother =
        decay_configuration_.particles_[single_decay->first];

    for (auto decay_products_info = single_decay->second.begin();
        decay_products_info != single_decay->second.end();
        ++decay_products_info) {

      std::vector<ParticleStateInfo> decay_daughters(
          generateParticleStateInfoList(
              decay_products_info->particle_indices_));

      DecayTree temp_tree(*decay_tree);
      DecayNode mother_node;
      mother_node.state_info_ =
          decay_configuration_.particles_[single_decay->first];
      mother_node.strength_and_phase_ =
          decay_products_info->decay_strength_info_and_phase_;
      temp_tree.createDecay(mother_node, decay_daughters);
      new_decay_trees.push_back(temp_tree);
    }
  }
  return new_decay_trees;
}

void DecayTreeFactory::removeFaultyDecayTrees(
    std::vector<DecayTree>& decay_tree_list) const {
  std::vector<DecayTree> good_decay_trees;
  good_decay_trees.reserve(decay_tree_list.size());

  for (unsigned int i = 0; i < decay_tree_list.size(); ++i) {
    if (!decay_tree_list[i].hasCycles()
        && !decay_tree_list[i].isDisconnected()) {
      good_decay_trees.push_back(decay_tree_list[i]);
    }
  }

  removeUndistinguishableCombinations(good_decay_trees);
  decay_tree_list = good_decay_trees;
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
        decay_tree.getHelicityDecayTree()[*vp.first].state_info_.id_information_.name_);
  }
  std::sort(name_list.begin(), name_list.end());
  return name_list;
}

} /* namespace HelicityFormalism */
