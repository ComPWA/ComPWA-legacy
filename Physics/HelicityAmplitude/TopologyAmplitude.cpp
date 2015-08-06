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

#include "TopologyAmplitude.hpp"

namespace HelicityFormalism {

TopologyAmplitude::TopologyAmplitude() {
  // TODO Auto-generated constructor stub

}

TopologyAmplitude::~TopologyAmplitude() {
  // TODO Auto-generated destructor stub
}

std::complex<double> TopologyAmplitude::evaluate(
    const std::vector<KinematicVariables>& kinematic_variables,
    const IndexList& evaluation_index_list) const {
  std::complex<double> result(0.0, 0.0);
  // loop over the list of concrete versions of sequential decays
  for (unsigned int sequential_decay_index = 0;
      sequential_decay_index < sequential_decay_amplitude_list_.size();
      ++sequential_decay_index) {
    std::complex<double> sequential_decay_result(1.0, 0.0);
    // loop over all the decay amplitudes within each sequential decay
    for (unsigned int two_body_decay_index = 0;
        two_body_decay_index
            < sequential_decay_amplitude_list_[sequential_decay_index].size();
        ++two_body_decay_index) {
      // the results for each amplitude evaluation are multiplied to the sequential decay result
      const KinematicVariables& single_decay_kinematic_variables =
          kinematic_variables[evaluation_index_list[two_body_decay_index]];
      sequential_decay_result *=
          sequential_decay_amplitude_list_[sequential_decay_index][two_body_decay_index].first->evaluate(
              single_decay_kinematic_variables)
              * sequential_decay_amplitude_list_[sequential_decay_index][two_body_decay_index].second->evaluate();
    }
    // the sequential decay results are just added
    result += sequential_decay_result;
  }
  return result;
}

} /* namespace HelicityFormalism */
