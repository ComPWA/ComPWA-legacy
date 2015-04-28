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

#include "HelicityAmplitude.hpp"

#include <cmath>

namespace HelicityFormalism {

HelicityAmplitude::HelicityAmplitude(const ParticleState& initial_state,
    const ParticleStatePair& final_state) :
    initial_state_(initial_state), final_state_(final_state) {
}

HelicityAmplitude::~HelicityAmplitude() {
}

void HelicityAmplitude::init() {
  spin_factor_ = sqrt((2 * initial_state_.J_ + 1) / (4 * M_PI));
}

std::complex<double> HelicityAmplitude::evaluate(
    const Vector4<double>& boosted_4vector) const {
  double theta(boosted_4vector.Theta());
  double phi(boosted_4vector.Phi());
  return spin_factor_
      * Wigner_D(phi, theta, -phi, initial_state_.J_, initial_state_.M_,
          final_state_.first.M_ - final_state_.second.M_)
      * dynamics_function_->evaluate();
}

} /* namespace HelicityFormalism */
