// Vector4 template class definition file. -*- C++ -*-
/* Copyright 2008 Mike Williams (mwill@jlab.org)
 *
 * This file is part of qft++.
 *
 * qft++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * qft++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with qft++.  If not, see <http://www.gnu.org/licenses/>.
 */
// Author: Mike Williams
#ifndef _Vector4_H
#define _Vector4_H
//_____________________________________________________________________________
/** @file Vector4.h
 *  @brief Vector4 template class definition file.
 */
//_____________________________________________________________________________
#include "Core/FourMomentum.hpp"
#include "Tensor.h"
#include <cmath>
namespace ComPWA {
namespace QFT {

//_____________________________________________________________________________
/** @class Vector4
 *  @brief \f$ x_{\mu}, p_{\mu} \f$ : 4-vectors and 4-momenta.
 *
 *  @author Mike Williams
 *
 * This class provides extra functionality needed for 4-vectors. The class
 * inherits from Tensor and in most cases is used as such.
 */
//_____________________________________________________________________________

template <typename _Tp> class Vector4 : public Tensor<_Tp> {

public:
  // create/copy/destroy:

  Vector4() : Tensor<_Tp>::Tensor(1) { /** Default Constructor */
  }

  /// Constructor (initialize the 4-vector to be (t,x,y,z))
  Vector4(typename Type<_Tp>::ParamType __t, typename Type<_Tp>::ParamType __x,
          typename Type<_Tp>::ParamType __y, typename Type<_Tp>::ParamType __z)
      : Tensor<_Tp>::Tensor(1) {
    this->SetV4(__t, __x, __y, __z);
  }

  Vector4(std::vector<_Tp> __v) : Tensor<_Tp>::Tensor(1) {
    if (__v.size() != 4)
      throw std::runtime_error(
          "Vector4::Vector4 | Vector doen't seem to be a 4-Vector.`");
    this->SetV4(__v.at(0), __v.at(1), __v.at(2), __v.at(3));
  }

  Vector4(FourMomentum __v) : Tensor<_Tp>::Tensor(1) {
    this->SetV4(__v.e(), __v.px(), __v.py(), __v.pz());
  }

  /// Copy Constructor
  Vector4(const Vector4<_Tp> &__v4) : Tensor<_Tp>::Tensor(__v4) {}

  virtual ~Vector4() { /** Destructor */
  }

  // Setters:

  /// Set the 4-vector to (t,x,y,z)
  void SetV4(typename Type<_Tp>::ParamType __t,
             typename Type<_Tp>::ParamType __x,
             typename Type<_Tp>::ParamType __y,
             typename Type<_Tp>::ParamType __z) {
    this->Element(0) = __t;
    this->Element(1) = __x;
    this->Element(2) = __y;
    this->Element(3) = __z;
  }

  /// Set the 4-vector to (e,px,py,pz)
  inline void SetP4(typename Type<_Tp>::ParamType __e,
                    typename Type<_Tp>::ParamType __px,
                    typename Type<_Tp>::ParamType __py,
                    typename Type<_Tp>::ParamType __pz) {
    this->SetV4(__e, __px, __py, __pz);
  }

  // Getters:

  /// Get the time component of the 4-vector
  inline const _Tp &T() const { return this->Element(0); }

  /// Get the x component of the 4-vector
  inline const _Tp &X() const { return this->Element(1); }

  /// Get the y component of the 4-vector
  inline const _Tp &Y() const { return this->Element(2); }

  /// Get the z component of the 4-vector
  inline const _Tp &Z() const { return this->Element(3); }

  /// Get the time component of the 4-vector
  inline _Tp &T() { return this->Element(0); }

  /// Get the x component of the 4-vector
  inline _Tp &X() { return this->Element(1); }

  /// Get the y component of the 4-vector
  inline _Tp &Y() { return this->Element(2); }

  /// Get the z component of the 4-vector
  inline _Tp &Z() { return this->Element(3); }

  /// Get the energy component of the 4-vector
  inline const _Tp &E() const { return this->Element(0); }

  /// Get the x component of the 4-vector
  inline const _Tp &Px() const { return this->Element(1); }

  /// Get the y component of the 4-vector
  inline const _Tp &Py() const { return this->Element(2); }

  /// Get the z component of the 4-vector
  inline const _Tp &Pz() const { return this->Element(3); }

  /// Get the energy component of the 4-vector
  inline _Tp &E() { return this->Element(0); }

  /// Get the x component of the 4-vector
  inline _Tp &Px() { return this->Element(1); }

  /// Get the y component of the 4-vector
  inline _Tp &Py() { return this->Element(2); }

  /// Get the z component of the 4-vector
  inline _Tp &Pz() { return this->Element(3); }

  // operators:

  /// Assignment operator
  template <typename T> Vector4<_Tp> &operator=(const Tensor<T> &__tensor) {
    if (__tensor.Rank() != 1) {
      std::cout
          << "Error! Attempt to set Vector4 equal to a tensor w/ rank != 1"
          << std::endl;
    }
    assert(__tensor.Rank() == 1);
    this->SetV4(__tensor(0), __tensor(1), __tensor(2), __tensor(3));

    return *this;
  }

  /// Sets @a this = @a this + @a v4
  template <typename T> Vector4<_Tp> &operator+=(const Vector4<T> &__v4) {
    *this = (*this) + __v4;
    return *this;
  }

  /// Sets @a this = @a this - @a v4
  template <typename T> Vector4<_Tp> &operator-=(const Vector4<T> &__v4) {
    *this = (*this) - __v4;
    return *this;
  }

  // Functions:

  /** Is this a valid 4-momentum?
   * Returns @a true if this is a valid 4-momentum. Valid means that all
   * of its elements are real numbers and its \f$ \beta \leq 1 \f$.
   */
  bool IsP4() const {
    if ((imag(this->E()) != 0.) || (imag(this->Px()) != 0.) ||
        (imag(this->Py()) != 0.) || (imag(this->Pz()) != 0.))
      return false;
    if (this->Beta() > 1.0)
      return false;
    return true;
  }

  /// Returns the magnitude of the momentum (\f$\sqrt{px^2 + py^2 + pz^2}\f$)
  inline _Tp P() const {
    return sqrt((this->Px()) * (this->Px()) + (this->Py()) * (this->Py()) +
                (this->Pz()) * (this->Pz()));
  }

  /// Returns \f$ \beta = \frac{p}{E} \f$
  inline _Tp Beta() const { return (this->P() / this->E()); }

  /// Returns \f$ \gamma = \frac{1}{\sqrt{1 - \beta^2}} \f$
  inline _Tp Gamma() const {
    return (1.0 / sqrt(1.0 - (this->Beta() * this->Beta())));
  }

  /// Returns \f$ \sqrt{p_{\mu} p^{\mu}} \f$
  inline _Tp Mass() const { return sqrt(this->Mass2()); }

  /// Returns \f$ p_{\mu} p^{\mu} \f$
  inline _Tp Mass2() const { return (*this) * (*this); }

  /// Returns \f$ \sqrt{p_{\mu} p^{\mu}} \f$
  inline _Tp M() const { return sqrt(this->Mass2()); }

  /// Returns \f$ p_{\mu} p^{\mu} \f$
  inline _Tp M2() const { return (*this) * (*this); }

  /// Returns the magnitude of the 4-vector (\f$\sqrt{x^2 + y^2 + z^2}\f$)
  inline _Tp R() const { return this->P(); }

  /// Returns \f$ \sqrt{x_{\mu} x^{\mu}} \f$
  inline _Tp Length() const { return this->Mass(); }

  /// Returns \f$ x_{\mu} x^{\mu} \f$
  inline _Tp Length2() const { return (*this) * (*this); }

  /// Returns \f$ \sqrt{x_{\mu} x^{\mu}} \f$
  inline _Tp L() const { return this->Mass(); }

  /// Returns \f$ x_{\mu} x^{\mu} \f$
  inline _Tp L2() const { return (*this) * (*this); }

  /// Returns \f$ cos(\theta) = \frac{z}{r} \f$
  inline _Tp CosTheta() const { return this->Z() / this->R(); }

  /// Returns \f$\sqrt{x^2 + y^2}\f$
  inline _Tp Rho() const {
    return sqrt((this->Px()) * (this->Px()) + (this->Py()) * (this->Py()));
  }

  /// Returns \f$\sqrt{px^2 + py^2}\f$
  inline _Tp Pxy() const {
    return sqrt((this->Px()) * (this->Px()) + (this->Py()) * (this->Py()));
  }

  /// Returns \f$ \theta = cos^{-1}(\frac{z}{r}) \f$
  inline _Tp Theta() const {
    if (this->Z() == 0)
      return acos(0); // return value is nan for Z=0 and R=0
    return acos(this->Z() / this->R());
  }

  /// Returns \f$ \theta = tan^{-1}(\frac{y}{x}) \f$  \f$ (-\pi,\pi) \f$
  inline _Tp Phi() const { return atan2(this->Y(), this->X()); }
};
//_____________________________________________________________________________
//
// some non-member operaors used to keep the type Vector4
//_____________________________________________________________________________

/// 4-vector addition
template <typename T1, typename T2>
Vector4<typename AddType<T1, T2>::Type> operator+(const Vector4<T1> &__v4a,
                                                  const Vector4<T2> &__v4b) {
  Vector4<typename AddType<T1, T2>::Type> ret;
  ret = __v4a.operator+(__v4b);
  return ret;
}
//_____________________________________________________________________________

/// 4-vector subtraction
template <typename T1, typename T2>
Vector4<typename SubType<T1, T2>::Type> operator-(const Vector4<T1> &__v4a,
                                                  const Vector4<T2> &__v4b) {
  Vector4<typename SubType<T1, T2>::Type> ret;
  ret = __v4a.operator-(__v4b);
  return ret;
}
//_____________________________________________________________________________

/// 4-vector contraction (returns the type)
template <typename T1, typename T2>
typename MultType<T1, T2>::Type operator*(const Vector4<T1> &__v4a,
                                          const Vector4<T2> &__v4b) {
  return (__v4a.operator*(__v4b)).Element();
}
//_____________________________________________________________________________

} // namespace QFT
} // namespace ComPWA
#endif /* _Vector4_H  */
