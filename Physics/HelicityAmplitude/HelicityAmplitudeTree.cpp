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
    const std::vector<HelicityAngles>& helicity_angles) const {
  std::complex<double> result(0.0, 0.0);
  // loop over the list of concrete versions of sequential decays
  for (unsigned int sequential_decay_index = 0;
      sequential_decay_index < amplitude_nodes_.size();
      ++sequential_decay_index) {
    std::complex<double> sequential_decay_result(1.0, 0.0);
    // loop over all the decay amplitudes within each sequential decay
    for (unsigned int i = 0;
        i < sequential_decay_amplitude_list_[sequential_decay_index].size();
        ++i) {
      // the results for each amplitude evaluation are multiplied to the sequential decay result
      sequential_decay_result *=
          sequential_decay_amplitude_list_[sequential_decay_index][i]->evaluate(
              helicity_angles[i]);
    }
    // the sequential decay results are just added
    result += sequential_decay_result;
  }
  return result;
}

/*std::complex<double> HelicityAmplitudeTree::evaluateIntermediateNode(
 unsigned int node_index,
 const std::vector<HelicityAngles>& helicity_angles) const {
 std::complex<double> result(0.0, 0.0);

 // if this node has links evaluate them first
 if (links_.size() > node_index) {
 // now evaluate them
 result = evaluateIntermediateNode(links_[node_index].first,
 helicity_angles);
 result *= evaluateIntermediateNode(links_[node_index].second,
 helicity_angles);
 }
 else {
 // then use that coordinate system to calculate the amplitude for this decay
 result = amplitude_nodes_[node_index]->evaluate(
 helicity_angles[node_index]);
 }

 return result;
 }*/

} /* namespace HelicityFormalism */
