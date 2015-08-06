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

#ifndef PHYSICS_HELICITYAMPLITUDE_TWOBODYDECAYAMPLITUDE_HPP_
#define PHYSICS_HELICITYAMPLITUDE_TWOBODYDECAYAMPLITUDE_HPP_

#include "ParticleStateDefinitions.hpp"
#include "AbstractDynamicalFunction.hpp"

#include <memory>
#include <complex>

namespace HelicityFormalism {

class TwoBodyDecayAmplitude {
  TwoBodyDecayInformation decay_info_;

  double spin_factor_;

public:
  TwoBodyDecayAmplitude(const TwoBodyDecayInformation& decay_info);
  virtual ~TwoBodyDecayAmplitude();

  void init();

  std::complex<double> evaluate(const KinematicVariables& helicity_angles) const;
};

} /* namespace HelicityFormalism */
#endif /* PHYSICS_HELICITYAMPLITUDE_TWOBODYDECAYAMPLITUDE_HPP_ */
