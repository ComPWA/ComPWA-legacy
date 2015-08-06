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

  bool isNodeADaughter(ParticleIndexDecayTree::const_iterator& node,
      const std::vector<unsigned int>& list_of_daughters) const;

  bool isNodeADaughterInTopology(ParticleIndexDecayTree::const_iterator& node,
      const ParticleIndexDecayTree& decay_topology) const;

  ParticleIndexDecayTree::const_iterator determineTopNode(
      const ParticleIndexDecayTree& decay_topology) const;

  bool canDecayTreesGrow(const std::vector<DecayTree>& decay_trees,
      const ParticleIndexDecayTree& decay_topology) const;

  std::vector<DecayTree> growSingleLeafOnDecayTrees(
      const std::vector<DecayTree>& decay_trees,
      ParticleIndexDecayTree::const_iterator& single_decay) const;

  std::vector<DecayTree> growNextDecayTreeLayer(
      const std::vector<DecayTree>& decay_trees,
      const ParticleIndexDecayTree& decay_topology) const;

  std::vector<DecayTree> createDecayTreeSeedList(
      ParticleIndexDecayTree::const_iterator& top_node_iter) const;

public:
  DecayTreeFactory(
      const DecayConfiguration& decay_configuration);
  virtual ~DecayTreeFactory();

  std::vector<DecayTree> createDecayTrees() const;
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_HELICITYDECAYTREEFACTORY_HPP_ */
