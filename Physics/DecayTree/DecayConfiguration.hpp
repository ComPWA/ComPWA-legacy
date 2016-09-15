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

#ifndef PHYSICS_HELICITYAMPLITUDE_DECAYCONFIGURATION_HPP_
#define PHYSICS_HELICITYAMPLITUDE_DECAYCONFIGURATION_HPP_

#include <map>
#include <vector>

#include "Physics/HelicityAmplitude/ParticleStateDefinitions.hpp"

namespace ComPWA {
namespace Physics {
namespace DecayTree {

struct DecayProductsInfo {
  std::vector<unsigned int> particle_indices_;
  boost::property_tree::ptree decay_strength_info_and_phase_;
};

typedef std::map<unsigned int, DecayProductsInfo> ParticleIndexDecayTree;

class DecayConfiguration {
  friend class DecayTreeFactory;

  std::vector<ParticleStateInfo> particles_;

  std::vector<ParticleIndexDecayTree> concrete_decay_trees_;

  ParticleIndexDecayTree current_concrete_decay_tree_;


  boost::property_tree::ptree createPropertyTreeForParticleIndexTree(
      const ParticleIndexDecayTree& particle_index_decay_tree) const;

  void addNextDecayTreeLayerToPropertyTree(
      boost::property_tree::ptree& decay_tree_pt,
      const std::vector<ParticleIndexDecayTree::const_iterator> & current_nodes) const;

  boost::property_tree::ptree createSingleDecayNodeToPropertyTree(
      const ParticleIndexDecayTree::const_iterator& current_node) const;

  void fillParticleLists(
      const ParticleIndexDecayTree& particle_index_decay_tree,
      std::map<unsigned int, unsigned int>& unique_final_state_particle_index_list,
      std::map<unsigned int, unsigned int>& unique_intermediate_state_particle_index_list) const;

  boost::property_tree::ptree createPropertyTreeForFinalStateParticle(
      const ParticleStateInfo& particle) const;

  boost::property_tree::ptree createPropertyTreeForIntermediateStateParticle(
      const ParticleStateInfo& particle) const;

public:
  DecayConfiguration();
  virtual ~DecayConfiguration();

  void addCurrentDecayTreeToList();

  void addDecayToCurrentDecayTree(const ParticleStateInfo& mother,
      const std::vector<ParticleStateInfo>& daughter_states,
      const boost::property_tree::ptree& decay_strength_info_and_phase);

  std::vector<unsigned int> addParticlesToList(
      const std::vector<ParticleStateInfo>& particle_list);

  unsigned int addParticleToList(ParticleStateInfo particle);

  void setRemainingParticleProperties(ParticleStateInfo& particle) const;

  boost::property_tree::ptree exportConfigurationToPropertyTree() const;

  ParticleIndexDecayTree::const_iterator determineTopNode(
      const ParticleIndexDecayTree& decay_topology) const;

  bool isNodeADaughterInTopology(ParticleIndexDecayTree::const_iterator& node,
      const ParticleIndexDecayTree& decay_topology) const;

  bool isNodeADaughter(ParticleIndexDecayTree::const_iterator& node,
      const std::vector<unsigned int>& list_of_daughters) const;

  void printInfo() const;
  void printDecayTree(const ParticleIndexDecayTree& index_decay_tree) const;

};

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_DECAYCONFIGURATION_HPP_ */
