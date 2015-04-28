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

#ifndef HELICITYSTATEDEFINITIONS_HPP_
#define HELICITYSTATEDEFINITIONS_HPP_

#include <qft++.h>

namespace HelicityFormalism {

struct Helicity4Vector {
  Vector4<double> total_momentum_;

  void boostIntoCMSOf(const Vector4<double>& mother_state) {
    PolVector state_vector;
    state_vector.SetP4(total_momentum_, total_momentum_.Mass());
    state_vector.Boost(mother_state);
    total_momentum_ = state_vector.GetP4();
  }
};

struct ParticleState {
  int particle_id_;
  Spin J_;
  Spin M_;
};

typedef std::pair<ParticleState, ParticleState> ParticleStatePair;
typedef std::pair<unsigned int, unsigned int> IndexPair;

}

#endif /* HELICITYSTATEDEFINITIONS_HPP_ */
