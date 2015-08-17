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
    const TwoBodyDecaySpinInformation& decay_info) :
    decay_info_(decay_info) {
  init();
}

TwoBodyDecayAmplitude::~TwoBodyDecayAmplitude() {
}

void TwoBodyDecayAmplitude::init() {
  spin_factor_ = sqrt((2 * decay_info_.initial_state_.J_ + 1) / (4 * M_PI));
}

std::complex<double> TwoBodyDecayAmplitude::evaluate(
    const HelicityAngles& helicity_angles) const {
  return spin_factor_
      * Wigner_D(helicity_angles.phi_, helicity_angles.theta_,
          -helicity_angles.phi_, decay_info_.initial_state_.J_,
          decay_info_.initial_state_.M_,
          decay_info_.final_state_.first.M_
              - decay_info_.final_state_.second.M_);
}

} /* namespace HelicityFormalism */
