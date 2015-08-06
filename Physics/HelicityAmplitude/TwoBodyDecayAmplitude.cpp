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

#include "TwoBodyDecayAmplitude.hpp"

#include <cmath>

namespace HelicityFormalism {

TwoBodyDecayAmplitude::TwoBodyDecayAmplitude(
    const TwoBodyDecayInformation& decay_info) :
    decay_info_(decay_info) {
  init();
}

TwoBodyDecayAmplitude::~TwoBodyDecayAmplitude() {
}

void TwoBodyDecayAmplitude::init() {
  spin_factor_ = sqrt((2 * decay_info_.initial_state_.J_ + 1) / (4 * M_PI));
}

std::complex<double> TwoBodyDecayAmplitude::evaluate(
    const KinematicVariables& kinematic_variables) const {
  return spin_factor_
      * Wigner_D(kinematic_variables.phi_, kinematic_variables.theta_,
          -kinematic_variables.phi_, decay_info_.initial_state_.J_,
          decay_info_.initial_state_.M_,
          decay_info_.final_state_.first.M_
              - decay_info_.final_state_.second.M_);
}

} /* namespace HelicityFormalism */
