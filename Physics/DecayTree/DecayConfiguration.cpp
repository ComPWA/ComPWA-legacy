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

#include <algorithm>
#include <stdexcept>

#include "DecayConfiguration.hpp"
#include "Core/PhysConst.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayTree {

DecayConfiguration::DecayConfiguration() {
}

DecayConfiguration::~DecayConfiguration() {
}

void DecayConfiguration::addCurrentDecayTreeToList() {
  concrete_decay_trees_.push_back(current_concrete_decay_tree_);
  current_concrete_decay_tree_.clear();
}

void DecayConfiguration::addDecayToCurrentDecayTree(
    const ParticleStateInfo& mother,
    const std::vector<ParticleStateInfo>& daughter_states,
    const boost::property_tree::ptree& decay_strength_info_and_phase) {

  // add particles to list if not already existent and get index
  unsigned int mother_state_index = addParticleToList(mother);

  DecayProductsInfo products;
  products.decay_strength_info_and_phase_ = decay_strength_info_and_phase;
  products.particle_indices_ = addParticlesToList(daughter_states);

  current_concrete_decay_tree_[mother_state_index] = products;
}

std::vector<unsigned int> DecayConfiguration::addParticlesToList(
    const std::vector<ParticleStateInfo>& particle_list) {
  std::vector<unsigned int> index_list;
  for (unsigned int i = 0; i < particle_list.size(); ++i) {
    index_list.push_back(addParticleToList(particle_list[i]));
  }
  return index_list;
}

unsigned int DecayConfiguration::addParticleToList(ParticleStateInfo particle) {
  // first make sure the contents of the particle are all set (correctly)
  setRemainingParticleProperties(particle);

  unsigned int index(0);
  std::vector<ParticleStateInfo>::iterator found_particle = std::find(
      particles_.begin(), particles_.end(), particle);
  if (found_particle != particles_.end()) {
    index = std::distance(particles_.begin(), found_particle);
  }
  else {
    index = particles_.size();
    particles_.push_back(particle);
  }
  return index;
}

void DecayConfiguration::setRemainingParticleProperties(
    ParticleStateInfo& particle) const {

  ComPWA::PhysConst &physics_constants = ComPWA::PhysConst::Instance();

  if (physics_constants.particleExists(particle.pid_information_.name_)) {
    if (particle.pid_information_.particle_id_
        != physics_constants.findParticle(particle.pid_information_.name_).id_) {
      particle.pid_information_.particle_id_ = physics_constants.findParticle(
          particle.pid_information_.name_).id_;
    }
    if (!particle.spin_information_.equalMagnitude(
        physics_constants.findParticle(particle.pid_information_.name_).getSpinLikeQuantumNumber(
            QuantumNumberIDs::SPIN))) {
      particle.spin_information_ = physics_constants.findParticle(
          particle.pid_information_.name_).getSpinLikeQuantumNumber(
          QuantumNumberIDs::SPIN);
    }
  }
}

boost::property_tree::ptree DecayConfiguration::exportConfigurationToPropertyTree() const {
// Create an empty property tree object

  boost::property_tree::ptree pt;

  std::map<unsigned int, unsigned int> unique_final_state_particle_index_list;
  std::map<unsigned int, unsigned int> unique_intermediate_state_particle_index_list;

  for (auto particle_index_decay_tree = concrete_decay_trees_.begin();
      particle_index_decay_tree != concrete_decay_trees_.end();
      ++particle_index_decay_tree) {
    boost::property_tree::ptree decay_tree_property_tree =
        createPropertyTreeForParticleIndexTree(*particle_index_decay_tree);

    pt.add_child("DecayTreeSetup.DecayTrees", decay_tree_property_tree);

    // fill the sets with the particles of the final state (only unique id matter here)

    fillParticleLists(*particle_index_decay_tree,
        unique_final_state_particle_index_list,
        unique_intermediate_state_particle_index_list);
  }

  for (auto particle_iter = unique_final_state_particle_index_list.begin();
      particle_iter != unique_final_state_particle_index_list.end();
      ++particle_iter) {
    pt.add_child("DecayTreeSetup.FinalState.Particle",
        createPropertyTreeForFinalStateParticle(
            particles_[particle_iter->second]));
  }

  for (auto particle_iter =
      unique_intermediate_state_particle_index_list.begin();
      particle_iter != unique_intermediate_state_particle_index_list.end();
      ++particle_iter) {
    pt.add_child("DecayTreeSetup.IntermediateStates.Particle",
        createPropertyTreeForIntermediateStateParticle(
            particles_[particle_iter->second]));
  }

  return pt;
}

boost::property_tree::ptree DecayConfiguration::createPropertyTreeForParticleIndexTree(
    const ParticleIndexDecayTree& particle_index_decay_tree) const {
  boost::property_tree::ptree decay_tree_pt;

  std::vector<ParticleIndexDecayTree::const_iterator> current_nodes;
  std::vector<ParticleIndexDecayTree::const_iterator> new_current_nodes;
  current_nodes.push_back(determineTopNode(particle_index_decay_tree));

  while (current_nodes.size() > 0) {
    addNextDecayTreeLayerToPropertyTree(decay_tree_pt, current_nodes);

    new_current_nodes.clear();

    for (auto current_node = current_nodes.begin();
        current_node != current_nodes.end(); ++current_node) {
      for (auto daughter_index =
          (*current_node)->second.particle_indices_.begin();
          daughter_index != (*current_node)->second.particle_indices_.end();
          ++daughter_index) {
        auto result = particle_index_decay_tree.find(*daughter_index);
        if (result != particle_index_decay_tree.end()) {
          new_current_nodes.push_back(result);
        }
      }
    }

    current_nodes = new_current_nodes;
  }

  return decay_tree_pt;
}

void DecayConfiguration::addNextDecayTreeLayerToPropertyTree(
    boost::property_tree::ptree& decay_tree_pt,
    const std::vector<ParticleIndexDecayTree::const_iterator> & current_nodes) const {
  for (auto node = current_nodes.begin(); node != current_nodes.end(); ++node) {
    decay_tree_pt.add_child("DecayTree.DecayNode",
        createSingleDecayNodeToPropertyTree(*node));
  }
}

boost::property_tree::ptree DecayConfiguration::createSingleDecayNodeToPropertyTree(
    const ParticleIndexDecayTree::const_iterator& current_node) const {
  boost::property_tree::ptree decay_node_pt;
  boost::property_tree::ptree mother_pt;
  mother_pt.put("id", particles_[current_node->first].unique_id_);
  mother_pt.put("J_z_numerator",
      particles_[current_node->first].spin_information_.J_z_numerator_);
  mother_pt.put("J_z_denominator",
      particles_[current_node->first].spin_information_.J_denominator_);

  decay_node_pt.add_child("StrengthPhase",
      current_node->second.decay_strength_info_and_phase_);

  decay_node_pt.add_child("Mother.ParticleState", mother_pt);

  for (auto daughter_index = current_node->second.particle_indices_.begin();
      daughter_index != current_node->second.particle_indices_.end();
      ++daughter_index) {
    boost::property_tree::ptree daughter_pt;
    daughter_pt.put("id", particles_[*daughter_index].unique_id_);
    daughter_pt.put("J_z_numerator",
        particles_[*daughter_index].spin_information_.J_z_numerator_);
    daughter_pt.put("J_z_denominator",
        particles_[*daughter_index].spin_information_.J_denominator_);
    decay_node_pt.add_child("Daughters.ParticleState", daughter_pt);
  }
  return decay_node_pt;
}

void DecayConfiguration::fillParticleLists(
    const ParticleIndexDecayTree& particle_index_decay_tree,
    std::map<unsigned int, unsigned int>& unique_final_state_particle_index_list,
    std::map<unsigned int, unsigned int>& unique_intermediate_state_particle_index_list) const {

  std::vector<ParticleIndexDecayTree::const_iterator> current_nodes;
  std::vector<ParticleIndexDecayTree::const_iterator> new_current_nodes;
  current_nodes.push_back(determineTopNode(particle_index_decay_tree));

  while (current_nodes.size() > 0) {
    for (auto decay_node = current_nodes.begin();
        decay_node != current_nodes.end(); ++decay_node) {
      if (unique_intermediate_state_particle_index_list.find(
          particles_[(*decay_node)->first].unique_id_)
          == unique_intermediate_state_particle_index_list.end()) {
        unique_intermediate_state_particle_index_list[particles_[(*decay_node)->first].unique_id_] =
            (*decay_node)->first;
      }
    }

    new_current_nodes.clear();

    for (auto current_node = current_nodes.begin();
        current_node != current_nodes.end(); ++current_node) {
      for (auto daughter_index =
          (*current_node)->second.particle_indices_.begin();
          daughter_index != (*current_node)->second.particle_indices_.end();
          ++daughter_index) {
        auto result = particle_index_decay_tree.find(*daughter_index);
        if (result != particle_index_decay_tree.end()) {
          new_current_nodes.push_back(result);
        }
        else {
          if (unique_final_state_particle_index_list.find(
              particles_[*daughter_index].unique_id_)
              == unique_final_state_particle_index_list.end()) {
            unique_final_state_particle_index_list[particles_[*daughter_index].unique_id_] =
                *daughter_index;
          }
        }
      }
    }

    current_nodes = new_current_nodes;
  }
}

boost::property_tree::ptree DecayConfiguration::createPropertyTreeForFinalStateParticle(
    const ParticleStateInfo& particle) const {
  boost::property_tree::ptree pt;
  pt.put("id", particle.unique_id_);
  pt.put("name", particle.pid_information_.name_);
  return pt;
}

boost::property_tree::ptree DecayConfiguration::createPropertyTreeForIntermediateStateParticle(
    const ParticleStateInfo& particle) const {
  boost::property_tree::ptree pt;
  pt.put("id", particle.unique_id_);
  pt.put("name", particle.pid_information_.name_);

//spin info
  boost::property_tree::ptree spin_pt;
  spin_pt.put("J_numerator", particle.spin_information_.J_numerator_);
  spin_pt.put("J_denominator", particle.spin_information_.J_denominator_);
  pt.add_child("SpinInfo", spin_pt);

//dynamical info
  pt.add_child("DynamicalInfo", particle.dynamical_information_);

  return pt;
}

ParticleIndexDecayTree::const_iterator DecayConfiguration::determineTopNode(
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

bool DecayConfiguration::isNodeADaughterInTopology(
    ParticleIndexDecayTree::const_iterator& node,
    const ParticleIndexDecayTree& decay_topology) const {
  bool is_daughter(false);
  ParticleIndexDecayTree::const_iterator other_decay_node;
// check if this index appears in any of the daughter lists if not this is our top node
  for (other_decay_node = decay_topology.begin();
      other_decay_node != decay_topology.end(); ++other_decay_node) {
    is_daughter = isNodeADaughter(node,
        other_decay_node->second.particle_indices_);
    if (is_daughter)
      break;
  }
  return is_daughter;
}

bool DecayConfiguration::isNodeADaughter(
    ParticleIndexDecayTree::const_iterator& node,
    const std::vector<unsigned int>& list_of_daughters) const {
  for (unsigned int daugther_index = 0;
      daugther_index != list_of_daughters.size(); ++daugther_index) {
    if (list_of_daughters[daugther_index] == node->first)
      return true;
  }
  return false;
}

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */
