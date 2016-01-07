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

#include <cmath>

#include "Core/DataPoint.hpp"
#include "Core/Kinematics.hpp"
#include "Physics/HelicityAmplitude/TwoBodyDecayAmplitude.hpp"

namespace HelicityFormalism {

TwoBodyDecayAmplitude::TwoBodyDecayAmplitude(
    const TwoBodyDecaySpinInformation& decay_info) :
    decay_info_(decay_info), spin_factor_(1.0), index_theta_helicity_angle_(0), index_phi_helicity_angle_(
        0) {
  init();
}

TwoBodyDecayAmplitude::~TwoBodyDecayAmplitude() {
}

void TwoBodyDecayAmplitude::init() {
  spin_factor_ = sqrt(
      (2.0 * decay_info_.initial_state_.J_numerator_
          / decay_info_.initial_state_.J_denominator_ + 1) / (4.0 * M_PI));

  index_theta_helicity_angle_ = Kinematics::instance()->getVariableIndex(
      "helicity_angle_theta");
  index_phi_helicity_angle_ = Kinematics::instance()->getVariableIndex(
      "helicity_angle_phi");

  J_.SetSpin(decay_info_.initial_state_.J_numerator_,
      decay_info_.initial_state_.J_denominator_);
  M_.SetSpin(decay_info_.initial_state_.J_z_numerator_,
      decay_info_.initial_state_.J_denominator_);
  d1_M_.SetSpin(decay_info_.final_state_.first.J_z_numerator_,
      decay_info_.final_state_.first.J_denominator_);
  d2_M_.SetSpin(decay_info_.final_state_.second.J_z_numerator_,
      decay_info_.final_state_.second.J_denominator_);
}

std::complex<double> TwoBodyDecayAmplitude::evaluate(const dataPoint& point,
    unsigned int evaluation_index) const {
  double theta(point.getVal(evaluation_index + index_theta_helicity_angle_));
  double phi(point.getVal(evaluation_index + index_phi_helicity_angle_));

  return spin_factor_ * Wigner_D(phi, theta, -phi, J_, M_, d1_M_ - d2_M_);
}

} /* namespace HelicityFormalism */
