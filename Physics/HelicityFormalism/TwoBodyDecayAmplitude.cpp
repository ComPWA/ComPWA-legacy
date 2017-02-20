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
#include "Core/DataPointStorage.hpp"
#include "Physics/HelicityAmplitude/TwoBodyDecayAmplitude.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

namespace ComPWA {
namespace Physics {
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
      (2.0 * decay_info_.initial_state_.GetNumerator()
          / decay_info_.initial_state_.GetDenominator() + 1) / (4.0 * M_PI));

  index_theta_helicity_angle_ =
		  Kinematics::instance()->FindVariable("helicity_angle_theta");
  index_phi_helicity_angle_ =
		  Kinematics::instance()->FindVariable("helicity_angle_phi");

  J_=decay_info_.initial_state_;
  M_=ComPWA::Spin(decay_info_.initial_state_.GetZNumerator(),
      decay_info_.initial_state_.GetDenominator());
  auto d1_M_=ComPWA::Spin(decay_info_.final_state_.first.GetZNumerator(),
      decay_info_.final_state_.first.GetDenominator());
  auto d2_M_=ComPWA::Spin(decay_info_.final_state_.second.GetZNumerator(),
      decay_info_.final_state_.second.GetDenominator());
  daughters_delta_M_ = d1_M_ - d2_M_;
}

std::complex<double> TwoBodyDecayAmplitude::evaluate(const dataPoint& point,
    unsigned int evaluation_index) const {
  double theta(point.getVal(evaluation_index + index_theta_helicity_angle_));
  double phi(point.getVal(evaluation_index + index_phi_helicity_angle_));

//  return spin_factor_ * Wigner_D(phi, theta, -phi, J_, M_, daughters_delta_M_);
  return spin_factor_ * ComPWA::Physics::AmplitudeSum::AmpWigner2::dynamicalFunction(
		  cos(phi), cos(theta), cos(-phi), J_, M_ , daughters_delta_M_
  );
}

std::complex<double> TwoBodyDecayAmplitude::evaluate(unsigned int storage_index,
    unsigned int data_index, unsigned int evaluation_index) const
{
  double theta(
      DataPointStorage::Instance().getDataList(storage_index,
          evaluation_index + index_theta_helicity_angle_)[data_index]);
  double phi(
      DataPointStorage::Instance().getDataList(storage_index,
          evaluation_index + index_phi_helicity_angle_)[data_index]);

  return spin_factor_ * ComPWA::Physics::AmplitudeSum::AmpWigner2::dynamicalFunction(
		  cos(phi), cos(theta), cos(-phi), J_, M_ , daughters_delta_M_
  );
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
