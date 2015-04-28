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

#include "HelicityDecayTree.hpp"

namespace HelicityFormalism {

HelicityDecayTree::HelicityDecayTree() {
  // TODO Auto-generated constructor stub

}

HelicityDecayTree::~HelicityDecayTree() {
  // TODO Auto-generated destructor stub
}

void HelicityDecayTree::createDecay(const ParticleState& mother,
    const ParticleStatePair& daughters) {

  // check if the mother particle exists as a leaf currently
  // or one of the daughters is the current top node
  // if this is not the case then only nothing could be inside the tree
  // otherwise throw exception..

  // then make the correct inserts into the vector and link appropriately
}

} /* namespace HelicityFormalism */
