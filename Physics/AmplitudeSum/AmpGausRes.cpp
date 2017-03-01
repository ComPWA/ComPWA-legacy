//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpGausRes.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

AmpGausRes::AmpGausRes(const char *name, unsigned int varIdA,
                       std::shared_ptr<DoubleParameter> mag,
                       std::shared_ptr<DoubleParameter> phase,
                       std::shared_ptr<DoubleParameter> mass,
                       std::shared_ptr<DoubleParameter> width,
                       std::string mother, std::string particleA,
                       std::string particleB, int nCalls, normStyle nS)
    : AmpAbsDynamicalFunction(
          name, varIdA, 0, mag, phase, mass, ComPWA::Spin(0.), ComPWA::Spin(0.),
          ComPWA::Spin(0.), +1, 0, mother, particleA, particleB,
          formFactorType::noFormFactor, nCalls, nS),
      _width(width) {}

AmpGausRes::~AmpGausRes() {}

std::complex<double> AmpGausRes::EvaluateAmp(const dataPoint &point) const {

  double m0 = _mass->GetValue();
  double width = _width->GetValue();
  double m = point.getVal(_subSys);

  std::complex<double> gaus(
      GetNormalization() * exp(-1 * (m - m0) * (m - m0) / width / width / 2.),
      0);

  return gaus;
}
std::complex<double> AmpGausRes::Evaluate(const dataPoint &point) const {
  return EvaluateAmp(point) * EvaluateWignerD(point);
}

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */
