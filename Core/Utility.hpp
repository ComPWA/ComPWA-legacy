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


//! Utility defines basic structures

#ifndef CORE_UTILITY_HPP_
#define CORE_UTILITY_HPP_

namespace ComPWA {

struct Spin {
  unsigned int J_numerator_;
  unsigned int J_denominator_;
  int J_z_numerator_;

  Spin() :
      J_numerator_(0), J_denominator_(1), J_z_numerator_(0) {
  }

  bool operator==(const Spin &rhs) const {
    if (this->J_numerator_ != rhs.J_numerator_)
      return false;
    if (this->J_denominator_ != rhs.J_denominator_)
      return false;
    if (this->J_z_numerator_ != rhs.J_z_numerator_)
      return false;

    return true;
  }
  bool operator!=(const Spin &rhs) const {
    return !(*this == rhs);
  }

  bool operator<(const Spin &rhs) const {
    if (this->J_numerator_ < rhs.J_numerator_)
      return true;
    else if (this->J_numerator_ > rhs.J_numerator_)
      return false;
    if (this->J_denominator_ < rhs.J_denominator_)
      return true;
    else if (this->J_denominator_ > rhs.J_denominator_)
      return false;
    if (this->J_z_numerator_ < rhs.J_z_numerator_)
      return true;

    return false;
  }
  bool operator>(const Spin &rhs) const {
    return (rhs < *this);
  }
};

}

#endif /* CORE_UTILITY_HPP_ */
