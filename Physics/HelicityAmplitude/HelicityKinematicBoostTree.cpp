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

/**
 * Calculates the boosted 4 vector of the first particle from the two body decay
 * in the helicity formalism.
 */
/*Vector4 determineBoostedKinematicVariables(
    std::pair<Vector4, Vector4> two_body_state, Vector4 mother) {

  // define particle 1 of the two body decay
  PolVector particle1_4vector;
  particle1_4vector.SetP4(two_body_state.first);
  // define the two body state
  PolVector two_body_4vector;
  two_body_4vector.SetP4(two_body_state.first + two_body_state.second);
  // boost particle1 into the rest frame of the two body state
  particle1_4vector.Boost(two_body_rest_frame.GetP4());
  // then boost the two body state into the rest frame of its mother
  two_body_4vector.Boost(mother);
  // now determine the theta and phi values of the boosted particle1 vector
  // with respect to the boosted two body state



}*/

const Vector4<double>& HelicityKinematicBoostTree::getBoosted4Vector(
    unsigned int node_index) const {
  return boosted_4vectors_[node_index];
}

} /* namespace HelicityFormalism */
