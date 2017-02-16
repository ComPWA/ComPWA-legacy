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
#include "Physics/HelicityAmplitude/CoherentIntensity.hpp"
#include "Physics/HelicityAmplitude/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

const double CoherentIntensity::Intensity(const dataPoint &point) {
  return IntensityNoEff(point) * eff_->Evaluate(point);
}

const double CoherentIntensity::IntensityNoEff(const dataPoint &point) {
  double result = 0.0;
  for (auto i : _seqDecayAmps)
    result += (*i)->Evaluate(point);
  return result;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
