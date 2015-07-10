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

struct ParticleState {
  int particle_id_;
  std::string name_;
  Spin J_;
  Spin M_;

  bool operator==(const ParticleState &rhs) const {
    if (this->particle_id_ != rhs.particle_id_)
      return false;
    if (this->name_ != rhs.name_)
      return false;
    if (this->J_.Numerator() != rhs.J_.Numerator())
      return false;
    if (this->J_.Denominator() != rhs.J_.Denominator())
      return false;
    if (this->M_.Numerator() != rhs.M_.Numerator())
      return false;
    if (this->M_.Denominator() != rhs.M_.Denominator())
      return false;

    return true;
  }
};

typedef std::pair<ParticleState, ParticleState> ParticleStatePair;
typedef std::pair<unsigned int, unsigned int> IndexPair;

}

#endif /* HELICITYSTATEDEFINITIONS_HPP_ */
