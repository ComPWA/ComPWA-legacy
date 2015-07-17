//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//--------------------------------------------------------------------------------

#include "HelicityKinematicBoostTree.hpp"

namespace HelicityFormalism {

HelicityKinematicBoostTree::HelicityKinematicBoostTree() {
  // TODO Auto-generated constructor stub

}

HelicityKinematicBoostTree::~HelicityKinematicBoostTree() {
  // TODO Auto-generated destructor stub
}

void HelicityKinematicBoostTree::constructBoostingTree(
    const HelicityDecayTree& decay_tree, const dataPoint& event_data) {

  // construct the boosting tree from the decay tree and and the data of this event
  // that means for each event this has to be computed over and over....
  // also we should use not only one decay tree but all of them, since the same
  // boosting appears across multiple of the trees

}



const Vector4<double>& HelicityKinematicBoostTree::getBoosted4Vector(
    unsigned int node_index) const {
  return boosted_4vectors_[node_index];
}

} /* namespace HelicityFormalism */
