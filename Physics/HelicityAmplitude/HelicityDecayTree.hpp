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

#ifndef HELICITYDECAYTREE_HPP_
#define HELICITYDECAYTREE_HPP_

#include "HelicityStateDefinitions.hpp"

namespace HelicityFormalism {

class HelicityDecayTree {
  std::vector<ParticleState> nodes_;
  std::vector<std::pair<unsigned int, IndexPair> > links_;

public:
  HelicityDecayTree();
  virtual ~HelicityDecayTree();

  void createDecay(const ParticleState& mother,
      const ParticleStatePair& daughters);
};

} /* namespace HelicityFormalism */

#endif /* HELICITYDECAYTREE_HPP_ */
