//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//----------------------------------------------------------------------------------

#include "Physics/HelicityFormalism/AbstractDynamicalFunction.hpp"
#include "Tools/Integration.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

//! Integral
double AbstractDynamicalFunction::Integral() const {
  auto intAlg = ComPWA::Tools::IntegralByQuadrature<AbstractDynamicalFunction>(
      *this, _limits);
  double integral = intAlg.Integral();
  LOG(trace) << "AbstractDynamicalFunction::Integral() | Integral is "
             << integral << ".";
  return std::sqrt(integral);
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
