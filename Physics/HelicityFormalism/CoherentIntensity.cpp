//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//--------------------------------------------------------------------------------
#include "Core/Efficiency.hpp"
#include "Physics/HelicityFormalism/CoherentIntensity.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

double CoherentIntensity::Intensity(const dataPoint &point) const {
  return IntensityNoEff(point) * _eff->evaluate(point);
}

double CoherentIntensity::IntensityNoEff(const dataPoint &point) const {
  std::complex<double> result(0.0,0.0);
  for (auto i : _seqDecays)
    result += i->Evaluate(point);
  return std::norm(result);
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
