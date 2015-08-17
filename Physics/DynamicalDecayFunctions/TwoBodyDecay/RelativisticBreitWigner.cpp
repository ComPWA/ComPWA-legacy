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

#include "Physics/DynamicalDecayFunctions/TwoBodyDecay/RelativisticBreitWigner.hpp"
#include "Physics/DynamicalDecayFunctions/Kinematics.hpp"

namespace DynamicalFunctions {

RelativisticBreitWigner::RelativisticBreitWigner(const Spin& J) :
    resonance_width_(new DoubleParameter("width")), resonance_mass_(
        new DoubleParameter("mass")), meson_radius_(
        new DoubleParameter("meson_radius")), J_(J) {

  parameter_list_.AddParameter(resonance_width_);
  parameter_list_.AddParameter(resonance_mass_);
  parameter_list_.AddParameter(meson_radius_);
}

RelativisticBreitWigner::~RelativisticBreitWigner() {
}

void RelativisticBreitWigner::initialiseParameters() {
}

std::complex<double> RelativisticBreitWigner::evaluate(const dataPoint& point,
    unsigned int evaluation_index) {

  double mSq = point.getVal(evaluation_index++);
  double ma = point.getVal(evaluation_index++);
  double mb = point.getVal(evaluation_index++);

  double mR = resonance_mass_->GetValue();
  double width = resonance_width_->GetValue();
  unsigned int J(J_.Numerator() / J_.Denominator());
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

}
