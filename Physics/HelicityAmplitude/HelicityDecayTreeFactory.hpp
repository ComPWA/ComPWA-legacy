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

#ifndef HELICITYDECAYTREEFACTORY_HPP_
#define HELICITYDECAYTREEFACTORY_HPP_

#include "HelicityDecayTree.hpp"
#include "HelicityDecayConfiguration.hpp"

namespace HelicityFormalism {

class HelicityDecayTreeFactory {
  const HelicityDecayConfiguration& decay_configuration_;

  bool isNodeADaughter(DecayTopology::const_iterator& node,
      const std::vector<unsigned int>& list_of_daughters) const;

  bool isNodeADaughterInTopology(DecayTopology::const_iterator& node,
      const DecayTopology& decay_topology) const;

  DecayTopology::const_iterator determineTopNode(
      const DecayTopology& decay_topology) const;

  bool canDecayTreesGrow(const std::vector<HelicityDecayTree>& decay_trees,
      const DecayTopology& decay_topology) const;

  std::vector<HelicityDecayTree> growSingleLeafOnDecayTrees(
      const std::vector<HelicityDecayTree>& decay_trees,
      DecayTopology::const_iterator& single_decay) const;

  std::vector<HelicityDecayTree> growNextDecayTreeLayer(
      const std::vector<HelicityDecayTree>& decay_trees,
      const DecayTopology& decay_topology) const;

  std::vector<HelicityDecayTree> createDecayTreeSeedList(
      DecayTopology::const_iterator& top_node_iter) const;

public:
  HelicityDecayTreeFactory(
      const HelicityDecayConfiguration& decay_configuration);
  virtual ~HelicityDecayTreeFactory();

  std::vector<std::vector<HelicityDecayTree> > createDecayTrees() const;
};

} /* namespace HelicityFormalism */

#endif /* HELICITYDECAYTREEFACTORY_HPP_ */
