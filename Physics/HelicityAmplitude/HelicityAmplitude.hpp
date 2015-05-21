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

#ifndef HELICITYAMPLITUDE_HPP_
#define HELICITYAMPLITUDE_HPP_

#include "HelicityStateDefinitions.hpp"
#include "AbstractDynamicalFunction.hpp"

#include <memory>
#include <complex>

namespace HelicityFormalism {

class HelicityAmplitude {
  ParticleState initial_state_;
  ParticleStatePair final_state_;

  std::shared_ptr<AbstractDynamicalFunction> dynamics_function_;

  double spin_factor_;

public:
  HelicityAmplitude(const ParticleState& initial_state,
      const ParticleStatePair& final_state);
  virtual ~HelicityAmplitude();

  void init();

  std::complex<double> evaluate(const Vector4<double>& boosted_4vector) const;
};

} /* namespace HelicityFormalism */
#endif /* HELICITYAMPLITUDE_HPP_ */
