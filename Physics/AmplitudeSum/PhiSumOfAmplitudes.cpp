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

#include <cmath>

#include "Physics/AmplitudeSum/PhiSumOfAmplitudes.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

PhiSumOfAmplitudes::PhiSumOfAmplitudes(const char *name) : _name(name) {}

PhiSumOfAmplitudes::PhiSumOfAmplitudes(const PhiSumOfAmplitudes &other,
                                       const char *name)
    : _name(name), _pdfList(other._pdfList), _intList(other._intList),
      _phaseList(other._phaseList), _angList(other._angList) {}
void PhiSumOfAmplitudes::addBW(std::shared_ptr<AmpAbsDynamicalFunction> theRes,
                               std::shared_ptr<DoubleParameter> r,
                               std::shared_ptr<DoubleParameter> phi,
                               std::shared_ptr<AmpWigner2> theAng) {
  _pdfList.push_back(theRes);
  _intList.push_back(r);
  _phaseList.push_back(phi);
  _angList.push_back(theAng);
}

void PhiSumOfAmplitudes::addBW(std::shared_ptr<AmpAbsDynamicalFunction> theRes,
                               std::shared_ptr<DoubleParameter> r,
                               std::shared_ptr<DoubleParameter> phi) {
  _pdfList.push_back(theRes);
  _intList.push_back(r);
  _phaseList.push_back(phi);
  _angList.push_back(std::shared_ptr<AmpWigner2>(new AmpWigner2(1, Spin(0.0))));
}

double PhiSumOfAmplitudes::evaluate() const {
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
  std::complex<double> res(0, 0);

  return fabs(3.1416 + atan2(res.imag(), res.real()));
}

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */
