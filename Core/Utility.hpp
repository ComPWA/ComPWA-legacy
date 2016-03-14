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

#include <string>
#include <map>
#include <vector>

#include "boost/property_tree/ptree.hpp"

namespace ComPWA {

typedef std::vector<unsigned int> IndexList;
typedef std::pair<unsigned int, unsigned int> IndexPair;
typedef std::map<unsigned int, unsigned int> IndexMapping;

typedef boost::property_tree::ptree DynamicalInfo;

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

struct IDInfo {
  int particle_id_;
  std::string name_;

  bool operator==(const IDInfo &rhs) const {
    /* if (this->id_ != rhs.id_)
     return false;*/
    if (this->particle_id_ != rhs.particle_id_)
      return false;
    if (this->name_ != rhs.name_)
      return false;

    return true;
  }
  bool operator!=(const IDInfo &rhs) const {
    return !(*this == rhs);
  }

  bool operator<(const IDInfo &rhs) const {
    /*if (this->id_ < rhs.id_)
     return true;
     else if (this->id_ > rhs.id_)
     return false;*/

    return lessThenIgnoringID(*this, rhs);
  }

  static bool lessThenIgnoringID(const IDInfo &lhs, const IDInfo &rhs) {
    if (lhs.particle_id_ < rhs.particle_id_)
      return true;
    else if (lhs.particle_id_ > rhs.particle_id_)
      return false;
    if (lhs.name_ < rhs.name_)
      return true;

    return false;
  }

  bool operator>(const IDInfo &rhs) const {
    return (rhs < *this);
  }
};

struct ParticleStateInfo {
  unsigned int unique_id_;
  IDInfo pid_information_;
  Spin spin_information_;
  DynamicalInfo dynamical_information_;

  bool operator==(const ParticleStateInfo &rhs) const {
    if (this->unique_id_ != rhs.unique_id_)
      return false;
    if (this->pid_information_ != rhs.pid_information_)
      return false;
    if (this->spin_information_ != rhs.spin_information_)
      return false;

    return true;
  }

  bool operator!=(const ParticleStateInfo &rhs) const {
    return !((*this) == rhs);
  }

  bool operator<(const ParticleStateInfo &rhs) const {
    if (this->unique_id_ < rhs.unique_id_)
      return true;
    else if (this->unique_id_ > rhs.unique_id_)
      return false;
    if (this->pid_information_ < rhs.pid_information_)
      return true;
    else if (this->pid_information_ > rhs.pid_information_)
      return false;
    if (this->spin_information_ < rhs.spin_information_)
      return true;

    return false;
  }
  bool operator>(const ParticleStateInfo &rhs) const {
    return (rhs < *this);
  }
};

} /* namespace ComPWA */

#endif /* CORE_UTILITY_HPP_ */
