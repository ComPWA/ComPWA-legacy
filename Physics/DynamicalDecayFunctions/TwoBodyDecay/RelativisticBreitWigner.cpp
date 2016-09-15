//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//     Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>

#include "boost/property_tree/ptree.hpp"

#include "Core/DataPointStorage.hpp"
#include "Core/Kinematics.hpp"
#include "Physics/DynamicalDecayFunctions/TwoBodyDecay/RelativisticBreitWigner.hpp"
#include "Physics/DynamicalDecayFunctions/Kinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace DynamicalFunctions {

RelativisticBreitWigner::RelativisticBreitWigner(const ParticleStateInfo& psi,
    const ExternalParameters& external_parameters) {
  J_ = psi.spin_information_;
  resonance_width_.reset(
      new DoubleParameter(psi.pid_information_.name_ + "_width"));
  resonance_mass_.reset(
      new DoubleParameter(psi.pid_information_.name_ + "_mass"));
  meson_radius_.reset(
      new DoubleParameter(psi.pid_information_.name_ + "_meson_radius"));

  parameter_list_.AddParameter(resonance_mass_);
  parameter_list_.AddParameter(resonance_width_);
  parameter_list_.AddParameter(meson_radius_);

  index_cms_mass_squared_ = ComPWA::Kinematics::instance()->getVariableIndex(
      "cms_mass_squared");

  initialiseParameters(psi.dynamical_information_, external_parameters);
}

RelativisticBreitWigner::~RelativisticBreitWigner() {
}

void RelativisticBreitWigner::initialiseParameters(
    const boost::property_tree::ptree& parameter_info,
    const ExternalParameters& external_parameters) {
  resonance_mass_->SetValue(
      parameter_info.get_child("mass").get<double>("value"));
  if (parameter_info.get_child("mass").get<double>("fix"))
    resonance_mass_->SetParameterFixed();
  else
    resonance_mass_->SetParameterFree();

  resonance_mass_->SetMinValue(
      parameter_info.get_child("mass").get<double>("min"));
  resonance_mass_->SetMaxValue(
      parameter_info.get_child("mass").get<double>("max"));
  //resonance_mass_->SetUseBounds(false);

  resonance_width_->SetValue(
      parameter_info.get_child("width").get<double>("value"));
  if (parameter_info.get_child("width").get<double>("fix"))
    resonance_width_->SetParameterFixed();
  else
    resonance_width_->SetParameterFree();
  resonance_width_->SetMinValue(
      parameter_info.get_child("width").get<double>("min"));
  resonance_width_->SetMaxValue(
      parameter_info.get_child("width").get<double>("max"));
  //resonance_width_->SetUseBounds(false);

  auto meson_radius_pt = parameter_info.get_child("mesonRadius");
  if (meson_radius_pt.get_optional<double>("value")) {
    meson_radius_->SetValue(meson_radius_pt.get<double>("value"));
    if (meson_radius_pt.get<double>("fix"))
      meson_radius_->SetParameterFixed();
    else
      meson_radius_->SetParameterFree();
    meson_radius_->SetMinValue(meson_radius_pt.get<double>("min"));
    meson_radius_->SetMaxValue(meson_radius_pt.get<double>("max"));
  }
  else {
    meson_radius_->SetValue(parameter_info.get<double>("mesonRadius"));
    meson_radius_->SetParameterFixed();
  }

  // try to extract daughter masses from external parameters
  daughter1_mass_ = external_parameters.parameters_.GetDoubleParameter(
      external_parameters.parameter_name_mapping_.at("daughter1_mass"));
  daughter2_mass_ = external_parameters.parameters_.GetDoubleParameter(
      external_parameters.parameter_name_mapping_.at("daughter2_mass"));
}

std::complex<double> RelativisticBreitWigner::evaluate(const dataPoint& point,
    unsigned int evaluation_index) const {

  double mSq = point.getVal(evaluation_index + index_cms_mass_squared_);

  return evaluate(mSq);
}

std::complex<double> RelativisticBreitWigner::evaluate(
    unsigned int storage_index, unsigned int data_index,
    unsigned int evaluation_index) const {

  auto const& data_vec = DataPointStorage::Instance().getDataList(storage_index,
      evaluation_index + index_cms_mass_squared_);

  double mSq = data_vec[data_index];

  return evaluate(mSq);
}

std::complex<double> RelativisticBreitWigner::evaluate(double mSq) const {
  double mR = resonance_mass_->GetValue();
  double width = resonance_width_->GetValue();
  double ma = daughter1_mass_->GetValue();
  double mb = daughter2_mass_->GetValue();
  unsigned int J(J_.J_numerator_ / J_.J_denominator_);
  double mesonRadius = meson_radius_->GetValue();

  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  double barrier = Kinematics::FormFactor(sqrtS, ma, mb, J, mesonRadius)
      / Kinematics::FormFactor(mR, ma, mb, J, mesonRadius);
  std::complex<double> qTerm = std::pow(
      (Kinematics::phspFactor(sqrtS, ma, mb)
          / Kinematics::phspFactor(mR, ma, mb)) * mR / sqrtS, (2 * J + 1));
  //Calculate coupling constant to final state
  std::complex<double> g_final = Kinematics::widthToCoupling(mSq, mR, width, ma,
      mb, J, mesonRadius);

  //Coupling constant from production reaction. In case of a particle decay the production
  //coupling doesn't depend in energy since the CM energy is in the (RC) system fixed to the
  //mass of the decaying particle
  double g_production = 1;

  std::complex<double> denom = std::complex<double>(mR * mR - mSq, 0)
      + (-1.0) * i * sqrtS * (width * qTerm * barrier);

  std::complex<double> result = g_final * g_production / denom;

  if (result.real() != result.real() || result.imag() != result.imag()) {
    std::cout << "RelativisticBreitWigner::evaluate() | " << barrier << " "
        << mR << " " << mSq << " " << ma << " " << mb << std::endl;
    return 0;
  }
  return result;
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
