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
#include <iostream>
#include "Core/Exceptions.hpp"

namespace ComPWA {

typedef std::vector<unsigned int> IndexList;
typedef std::pair<unsigned int, unsigned int> IndexPair;
typedef std::map<unsigned int, unsigned int> IndexMapping;



class Spin {

public:
  double Val() const {
	  if (J_denominator_ == 0) return 0;
	  return (double)J_numerator_/J_denominator_;
  }

  Spin() :
      J_numerator_(0), J_denominator_(1), J_z_numerator_(0),
	  z_component_relevant(true) { }

  //! Constructor for half-integer spin
  Spin(int num, int denom) :
      J_numerator_(num), J_denominator_(denom),
	  J_z_numerator_(0), z_component_relevant(false) { }

  //! Constructor for integer spin
  Spin(int intSpin) :
      J_numerator_(intSpin), J_denominator_(1), J_z_numerator_(0),
	  z_component_relevant(true) { }

  bool equalMagnitude(const Spin &rhs) const {
    if (1.0 * this->J_numerator_ / this->J_denominator_
        == 1.0 * rhs.J_numerator_ / rhs.J_denominator_)
      return true;
    return false;
  }

  Spin operator+(const Spin &rhs) const {
	  /* We calculate (a/b - c/d) = (a*d - c*b)/(bd) */
	  int num = GetNumerator() * rhs.GetDenominator()
			  + rhs.GetNumerator() * GetDenominator();
	  int denom = rhs.GetDenominator() * GetDenominator();
	  int znum = GetZNumerator() * rhs.GetDenominator()
			  + rhs.GetZNumerator() * GetDenominator();
	  Spin s;
	  s.SetNumerator(num);
	  s.SetZNumerator(znum);
	  s.SetDenominator(denom);
	  s.Simplify();
	  return s;
  }

  Spin operator-(const Spin &rhs) const {
	  /* We calculate (a/b - c/d) = (a*d - c*b)/(bd) */
	  int num = GetNumerator() * rhs.GetDenominator()
			  - rhs.GetNumerator() * GetDenominator();
	  int denom = rhs.GetDenominator() * GetDenominator();
	  int znum = GetZNumerator() * rhs.GetDenominator()
			  - rhs.GetZNumerator() * GetDenominator();
	  Spin s;
	  s.SetNumerator(num);
	  s.SetZNumerator(znum);
	  s.SetDenominator(denom);
	  s.Simplify();
	  return s;
  }

  Spin& operator--()
  {
	  //Only decrement if the result is larger zero
	  if(GetNumerator()-1 >= 0){
		  SetNumerator( GetNumerator()-1 );
		  if(z_component_relevant)
			  SetZNumerator( GetZNumerator()+1 );
	  }
      return *this;
  }

  Spin operator--(int)
  {
      Spin tmp(*this); // copy
      operator++(); // pre-increment
      return tmp;   // return old value
  }

  Spin& operator++()
  {
	  SetNumerator( GetNumerator()+1 );
	  if(z_component_relevant)
		  SetZNumerator( GetZNumerator()+1 );
      return *this;
  }

  Spin operator++(int)
  {
      Spin tmp(*this); // copy
      operator++(); // pre-increment
      return tmp;   // return old value
  }

  void Simplify() {
	  int tmp, tmpZ;
	  tmp = tmpZ = ggT(GetNumerator(),GetDenominator());
	  if( UseZ() )
		  tmpZ = ggT(GetZNumerator(),GetDenominator());
	  if(tmp == tmpZ){
		  SetNumerator( GetNumerator()/tmp );
		  SetZNumerator( GetZNumerator()/tmp );
		  SetDenominator( GetDenominator()/tmp );
	  }
  }

  /**! Calculate largest common factor */
  static unsigned int ggT(unsigned int a, unsigned int b){
	  if(b == 0)
		  return a;
	  else return ggT(b, a % b);
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

  unsigned int GetNumerator() const { return J_numerator_; }

  unsigned int GetZNumerator() const { return J_z_numerator_; }

  unsigned int GetDenominator() const { return J_denominator_; }

  void SetNumerator(unsigned int num) {
	  J_numerator_ = num;
  }

  void SetZNumerator(unsigned int znum) {
	  J_z_numerator_ = znum;
  }

  void SetDenominator(unsigned int denom) {
	  if( denom != 1 && denom != 2)
		  throw BadParameter("Spin::SetDenominator() |"
				  " Should be equal 1 oder 2!");
	  J_denominator_ = denom;
  }

  double GetSpin() { return (double)J_numerator_/J_denominator_; }

  double GetZComponent() { return (double)J_z_numerator_/J_denominator_; }

  bool UseZ() const { return z_component_relevant; }

  void SetUseZ( bool b ) { z_component_relevant = b; }

protected:

  unsigned int J_numerator_;
  unsigned int J_denominator_;
  int J_z_numerator_;

  bool z_component_relevant;
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
  bool coherent;

  bool operator==(const ParticleStateInfo &rhs) const {
    if (this->unique_id_ != rhs.unique_id_)
      return false;
    if (this->pid_information_ != rhs.pid_information_)
      return false;
    if (this->spin_information_ != rhs.spin_information_)
      return false;
    if (this->coherent != rhs.coherent)
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
    if (this->coherent < rhs.coherent)
      return true;
    else if (this->coherent > rhs.coherent)
      return false;
    if (this->spin_information_ < rhs.spin_information_)
      return true;

    return false;
  }
  bool operator>(const ParticleStateInfo &rhs) const {
    return (rhs < *this);
  }

  friend std::ostream& operator<<(std::ostream& os,
      const ParticleStateInfo &rhs) {
    os << "unique id: " << rhs.unique_id_ << std::endl;
    os << "name: " << rhs.pid_information_.name_ << std::endl;
    os << "pid: " << rhs.pid_information_.particle_id_ << std::endl;
    os << "J: " << rhs.spin_information_.GetNumerator() << "/"
        << rhs.spin_information_.GetDenominator()<< "("
        << rhs.spin_information_.GetZNumerator()<< ")";
    if (rhs.coherent)
      os << " coherent" << std::endl;
    return os;
  }
};

} /* namespace ComPWA */

#endif /* CORE_UTILITY_HPP_ */
