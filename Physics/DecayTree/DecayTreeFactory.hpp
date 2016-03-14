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

namespace ComPWA {
namespace Physics {
namespace DecayTree {

class DecayTreeFactory {
  const DecayConfiguration& decay_configuration_;

  DecayTree createSeedDecayTree(
      ParticleIndexDecayTree::const_iterator & top_node_iter) const;

  void growSingleLeafOnDecayTree(DecayTree& decay_tree,
      ParticleIndexDecayTree::const_iterator& single_decay) const;

  std::vector<ParticleStateInfo> generateParticleStateInfoList(
      const std::vector<unsigned int>& particle_index_list) const;

  bool canDecayTreeGrow(const DecayTree& decay_tree,
      const ParticleIndexDecayTree& decay_topology) const;

  void growNextDecayTreeLayer(DecayTree& decay_tree,
      const ParticleIndexDecayTree& decay_topology) const;

  bool isDecayTreeFaulty(const DecayTree& decay_tree) const;

  void removeUndistinguishableCombinations(
      std::vector<DecayTree>& decay_tree_list) const;

  std::vector<std::string> convertToUniqueNameList(
      const DecayTree& decay_tree) const;

public:
  DecayTreeFactory(const DecayConfiguration& decay_configuration);
  virtual ~DecayTreeFactory();

  std::vector<DecayTree> createDecayTrees() const;
};

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_HELICITYDECAYTREEFACTORY_HPP_ */
