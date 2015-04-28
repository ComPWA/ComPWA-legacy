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

#include "HelicityAmplitudeTree.hpp"

namespace HelicityFormalism {

HelicityAmplitudeTree::HelicityAmplitudeTree() {
  // TODO Auto-generated constructor stub

}

HelicityAmplitudeTree::~HelicityAmplitudeTree() {
  // TODO Auto-generated destructor stub
}

std::complex<double> HelicityAmplitudeTree::evaluate(
    const HelicityKinematicBoostTree& boosted_4vectors) const {
  std::complex<double> result(0.0, 0.0);
  if (amplitude_nodes_.size() > 0) {
    // look for links of that node
    result = evaluateIntermediateNode(0, boosted_4vectors);
  }
  return result;
}

std::complex<double> HelicityAmplitudeTree::evaluateIntermediateNode(
    unsigned int node_index,
    const HelicityKinematicBoostTree& boosted_4vectors) const {
  std::complex<double> result(0.0, 0.0);

  // if this node has links evaluate them first
  if (links_.size() > node_index) {
    // now evaluate them
    result = evaluateIntermediateNode(links_[node_index].first,
        boosted_4vectors);
    result *= evaluateIntermediateNode(links_[node_index].second,
        boosted_4vectors);
  }
  else {
    // then use that coordinate system to calculate the amplitude for this decay
    result = amplitude_nodes_[node_index]->evaluate(
        boosted_4vectors.getBoosted4Vector(node_index));
  }

  return result;
}

} /* namespace HelicityFormalism */
