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

#include <complex>

#include "qft++.h"

#include "Physics/HelicityAmplitude/ParticleStateDefinitions.hpp"

namespace ComPWA {

class dataPoint;

namespace Physics {
namespace HelicityFormalism {

class TwoBodyDecayAmplitude {
  TwoBodyDecaySpinInformation decay_info_;

  double spin_factor_;

  unsigned int index_theta_helicity_angle_;
  unsigned int index_phi_helicity_angle_;

  ::Spin J_; // TODO: switch this to our definition eventually
  ::Spin M_;
  ::Spin d1_M_;
  ::Spin d2_M_;
  ::Spin daughters_delta_M_;

public:
  TwoBodyDecayAmplitude(const TwoBodyDecaySpinInformation& decay_info);
  virtual ~TwoBodyDecayAmplitude();

  void init();

  std::complex<double> evaluate(const dataPoint& point,
      unsigned int evaluation_index) const;

  std::complex<double> evaluate(unsigned int storage_index, unsigned int data_index,
      unsigned int evaluation_index) const;
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_HELICITYAMPLITUDE_TWOBODYDECAYAMPLITUDE_HPP_ */
