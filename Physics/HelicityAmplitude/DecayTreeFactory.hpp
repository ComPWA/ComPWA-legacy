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

#ifndef PHYSICS_HELICITYAMPLITUDE_DECAYTREEFACTORY_HPP_
#define PHYSICS_HELICITYAMPLITUDE_DECAYTREEFACTORY_HPP_

#include "DecayTree.hpp"
#include "DecayConfiguration.hpp"

namespace HelicityFormalism {

class DecayTreeFactory {
  const DecayConfiguration& decay_configuration_;

  ParticleIndexDecayTree::const_iterator determineTopNode(
      const ParticleIndexDecayTree& decay_topology) const;

  bool isNodeADaughterInTopology(ParticleIndexDecayTree::const_iterator& node,
      const ParticleIndexDecayTree& decay_topology) const;

  bool isNodeADaughter(ParticleIndexDecayTree::const_iterator& node,
      const std::vector<unsigned int>& list_of_daughters) const;

  std::vector<DecayTree> createDecayTreeSeedList(
      ParticleIndexDecayTree::const_iterator& top_node_iter) const;

  std::vector<std::vector<ParticleState> > makeDaughterCombinations(
      std::vector<std::vector<unsigned int> > particle_index_lists) const;

  void extendDaughterCombinations(
      std::vector<std::vector<unsigned int> > remaining_particle_index_lists,
      const std::vector<ParticleState>& decay_daughters_combination,
      std::vector<std::vector<ParticleState> >& decay_daughters_combination_list) const;

  bool isParticleValidForCombination(const ParticleState& particle,
      const std::vector<ParticleState>& combination) const;

  void removeUndistinguishableCombinations(
      std::vector<std::vector<ParticleState> >& decay_daughters_combination_lists) const;

  std::vector<unsigned int> convertToPIDList(
      const std::vector<ParticleState>& particle_list) const;

  bool canDecayTreesGrow(const std::vector<DecayTree>& decay_trees,
      const ParticleIndexDecayTree& decay_topology) const;

  std::vector<DecayTree> growNextDecayTreeLayer(
      const std::vector<DecayTree>& decay_trees,
      const ParticleIndexDecayTree& decay_topology) const;

  std::vector<DecayTree> growSingleLeafOnDecayTrees(
      const std::vector<DecayTree>& decay_trees,
      ParticleIndexDecayTree::const_iterator& single_decay) const;

public:
  DecayTreeFactory(const DecayConfiguration& decay_configuration);
  virtual ~DecayTreeFactory();

  std::vector<DecayTree> createDecayTrees() const;
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_HELICITYDECAYTREEFACTORY_HPP_ */
