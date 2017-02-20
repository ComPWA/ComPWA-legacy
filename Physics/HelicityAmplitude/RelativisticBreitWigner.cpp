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
#include "Physics/DynamicalDecayFunctions/RelativisticBreitWigner.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

std::complex<double> Evaluate(const dataPoint &point, int pos) const {
  std::complex<double> result = dynamicalFunction(
      point.getVal(pos), _mass->GetValue(), _mass1, _mass2, _width->GetValue(),
      (double)_spin, _mesonRadius->GetValue(), _ffType);
  return result;
}

std::complex<double> RelativisticBreitWigner::dynamicalFunction(
    double mSq, double mR, double ma, double mb, double width, unsigned int J,
    double mesonRadius, formFactorType ffType) {
  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  double barrier =
      Kinematics::FormFactor(sqrtS, ma, mb, J, mesonRadius, ffType) /
      Kinematics::FormFactor(mR, ma, mb, J, mesonRadius, ffType);

  std::complex<double> qTerm = std::pow((Kinematics::phspFactor(sqrtS, ma, mb) /
                                         Kinematics::phspFactor(mR, ma, mb)) *
                                            mR / sqrtS,
                                        (2 * J + 1));

  // Calculate coupling constant to final state
  std::complex<double> g_final =
      widthToCoupling(mSq, mR, width, ma, mb, J, mesonRadius, ffType);

  /*Coupling constant from production reaction. In case of a particle decay
   * the production coupling doesn't depend in energy since the CM energy
   * is in the (RC) system fixed to the mass of the decaying particle */
  double g_production = 1;

  std::complex<double> denom = std::complex<double>(mR * mR - mSq, 0) +
                               (-1.0) * i * sqrtS * (width * qTerm * barrier);

  std::complex<double> result = g_final * g_production / denom;

#ifndef NDEBUG
  if (std::isnan(result.real()) || std::isnan(result.imag())) {
    std::cout << "RelativisticBreitWigner::dynamicalFunction() | " << barrier
              << " " << mR << " " << mSq << " " << ma << " " << mb << std::endl;
    return 0;
  }
#endif

  return result;
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
