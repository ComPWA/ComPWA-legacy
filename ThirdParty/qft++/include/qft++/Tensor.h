// Tensor class definition file -*- C++ -*-
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
#ifndef _Tensor_H
#define _Tensor_H
//_____________________________________________________________________________
// Standard C++ Headers:
#include <iostream>
#include <cassert>
#include <vector>
#include <complex>
#include <cstdlib>
// Local Headers
#include "qft++/Type.h"
#include "qft++/OperationType.h"
#include "qft++/Conversion.h"
#include "qft++/SelectiveInclusion.h"
#include "qft++/Matrix_Base.h"
#include "qft++/Tensor_Base.h"
#include "qft++/TensorIndex.h"
//_____________________________________________________________________________
/** @file Tensor.h
 *  @brief Tensor template class definition file.
 */
//_____________________________________________________________________________
namespace ComPWA {
namespace QFT {

// Some utility functions templates:

/// Return @a zero for the type (defaults to numeric 0)
template <typename _Tp> inline _Tp zero(const _Tp &__var) { return 0; }

///// Return the conjugate for the type (defualts to the variable itself)
// template <typename _Tp> inline _Tp conj(const _Tp &__var){return __var;}

/// Return the conjugate for the type (defualts to the variable itself)
inline int conj(const int &__var) { return __var; }

/// Return the conjugate for the type (defualts to the variable itself)
inline double conj(const double &__var) { return __var; }

/// Return the imaginary part of the type (defaults to 0)
template <typename _Tp> inline _Tp imag(const _Tp &__var) { return 0; }

/// Return @a unity for the type (defaults to numeric 1)
template <typename _Tp> inline _Tp unity(const _Tp &__var) { return 1; }

//_____________________________________________________________________________
/** @class Tensor
 *  @author Mike Williams
 *
 *  @brief General template class for handling tensors and tensor operations.
 *
 * Tensor is a template class for handling tensors and tensor operations.
 * A Tensor object can be any rank and store any type that can be
 * stored in a C++ container class. This class has been written to be as
 * general and flexible as possible. Parameter passing has been optimized
 * using the Type class. Tensor has been designed so that two instantiations,
 * Tensor<_A> and Tensor<_B>, are completely compatible as long as _A and _B
 * are compatible (eg. _A*_B,_A+_B,...are defined).  Return types of
 * <em>mixed type</em> tensor operations are determined by the OperationType
 * template classes (MultType,DivType,AddType,SubType).
 *
 * <b>Example Usage </b>
 *
 * \include Tensor.ex
 */
//_____________________________________________________________________________

template <typename _Tp> class Tensor : public Tensor_Base {

protected:
  // attributes:
  std::vector<_Tp> _data; ///< Tensor elements (type _Tp)

private:
  // private functions:

  /// Copy @a tensor elements to @a this tensor
  template <typename T> void _Copy(const Tensor<T> &__tensor) {
    int size = __tensor.Size();
    _data.resize(size);
    for (int i = 0; i < size; i++)
      _data[i] = __tensor._data[i];
  }

protected:
  // friends:
  template <typename T> friend class Tensor;

public:
  // create/copy/destroy:

  /** Default Constructor (rank 0)*/
  Tensor() : Tensor_Base(), _data(1) {}

  /// Constructor
  /** @param rank Rank of the Tensor */
  Tensor(int __rank) : Tensor_Base(__rank), _data(1 << (__rank << 1)) {}

  /** Constructor
   * @param rank Rank of the Tensor
   *  @param init  Initial value of Tensor elements
   */
  Tensor(int __rank, typename Type<_Tp>::ParamType __init)
      : Tensor_Base(__rank), _data(1 << (__rank << 1), __init) {}

  /// Copy Constructor
  template <typename T>
  Tensor(const Tensor<T> &__tensor) : Tensor_Base(__tensor) {
    this->_Copy(__tensor);
  }

  /** Destructor */
  virtual ~Tensor() {}

  // basic functions:

  /// Returns the number of elements in the tensor
  inline int Size() const { return (1 << (_rank << 1)); }

  /** Set each element of @a this tensor to zero
   * Note: Legal if zero(_Tp) is legal
   *  (see c++-template-utils/TemplateUtilFuncs.h)
   */
  inline void Zero() {
    _Tp var_type;
    for (int i = 0; i < this->Size(); i++)
      _data[i] = zero(var_type);
  }

  /// Removes all elements from the tensor
  void Clear() {
    if (!_data.empty())
      _data.clear();
    _rank = -1;
  }

  /// Set the rank of the tensor to @a rank
  inline void SetRank(int __rank) {
    _data.resize(1 << (__rank << 1));
    _rank = __rank;
  }

  /** Boost using transformation defined by \f$\vec{\beta}=(bx,by,bz)\f$
   *
   *  Set \f$ X_{\mu\nu\ldots} =  X_{\delta\pi\ldots}
   *  \Lambda^{\delta}{}_{\mu}(\vec{\beta})\Lambda^{\pi}{}_{\nu}(\vec{\beta})
   *  \ldots \f$.
   */
  void Boost(double __bx, double __by, double __bz) {

    // check to see if bx,by,bz are all less than 1
    if (std::abs(__bx) >= 1 || std::abs(__by) >= 1 || std::abs(__bz) >= 1)
      std::cout << "Error! Attempt to boost using invalid boost vector."
                << std::endl;
    assert((std::abs(__bx) < 1) && (std::abs(__by) < 1) &&
           (std::abs(__bz) < 1));

    Tensor<double> lt(2); // Lorentz transformation tensor
    double gamma = 1.0 / sqrt(1.0 - __bx * __bx - __by * __by - __bz * __bz);
    double gamFact = (gamma * gamma) / (gamma + 1.0);

    // set up the Lorentz transformation tensor
    lt.Zero();

    lt(0, 0) = gamma;
    lt(0, 1) = gamma * __bx;
    lt(0, 2) = gamma * __by;
    lt(0, 3) = gamma * __bz;

    lt(1, 1) = (__bx * __bx * gamFact) + 1;
    lt(1, 2) = __bx * __by * gamFact;
    lt(1, 3) = __bx * __bz * gamFact;

    lt(2, 2) = (__by * __by * gamFact) + 1;
    lt(2, 3) = __by * __bz * gamFact;

    lt(3, 3) = (__bz * __bz * gamFact) + 1;

    lt(1, 0) = lt(0, 1);
    lt(2, 0) = lt(0, 2);
    lt(2, 1) = lt(1, 2);
    lt(3, 0) = lt(0, 3);
    lt(3, 1) = lt(1, 3);
    lt(3, 2) = lt(2, 3);

    this->Transform(lt);
  }

  /// Boost the Tensor to the rest frame of the 4-momentum @a p4.
  void Boost(const Tensor<double> &__p4) {
    if (__p4.Rank() != 1)
      std::cout << "Error! 4-momentum NOT rank 1." << std::endl;
    assert(__p4.Rank() == 1);

    this->Boost(-(__p4(1) / __p4(0)), -(__p4(2) / __p4(0)),
                -(__p4(3) / __p4(0)));
  }

  /// Rotate the tensor using Euler angles \f$ \alpha,\beta,\gamma \f$.
  void Rotate(double __alpha, double __beta, double __gamma) {

    double ca = cos(__alpha);
    double sa = sin(__alpha);
    double cb = cos(__beta);
    double sb = sin(__beta);
    double cg = cos(__gamma);
    double sg = sin(__gamma);

    Tensor<double> lt(2); // Lorentz transformation tensor

    lt.Zero();

    lt(0, 0) = 1.0;

    lt(1, 1) = ca * cb * cg - sa * sg;
    lt(1, 2) = sa * cb * cg + ca * sg;
    lt(1, 3) = -sb * cg;

    lt(2, 1) = -ca * cb * sg - sa * cg;
    lt(2, 2) = -sa * cb * sg + ca * cg;
    lt(2, 3) = sb * sg;

    lt(3, 1) = ca * sb;
    lt(3, 2) = sa * sb;
    lt(3, 3) = cb;

    this->Transform(lt);
  }

  /// Rotate about the x-axis
  void RotateX(double __alpha) {
    double ca = cos(__alpha);
    double sa = sin(__alpha);
    Tensor<double> lt(2); // Lorentz transformation tensor
    lt.Zero();

    lt(0, 0) = 1.0;
    lt(1, 1) = 1.0;
    lt(2, 2) = ca;
    lt(2, 3) = -sa;
    lt(3, 2) = sa;
    lt(3, 3) = ca;

    this->Transform(lt);
  }
  //_____________________________________________________________________________

  /// Rotate about the y-axis
  void RotateY(double __alpha) {

    double ca = cos(__alpha);
    double sa = sin(__alpha);
    Tensor<double> lt(2); // Lorentz transformation tensor
    lt.Zero();

    lt(0, 0) = 1.0;
    lt(1, 1) = ca;
    lt(1, 3) = sa;
    lt(2, 2) = 1.0;
    lt(3, 1) = -sa;
    lt(3, 3) = ca;

    this->Transform(lt);
  }
  //_____________________________________________________________________________

  /// Rotate about the z-axis
  void RotateZ(double __alpha) {
    double ca = cos(__alpha);
    double sa = sin(__alpha);
    Tensor<double> lt(2); // Lorentz transformation tensor
    lt.Zero();

    lt(0, 0) = 1.0;
    lt(1, 1) = ca;
    lt(1, 2) = -sa;
    lt(2, 1) = sa;
    lt(2, 2) = ca;
    lt(3, 3) = 1.0;

    this->Transform(lt);
  }

  /** Send the values of the tensor elements to @a os.
   *
   *  @param os ostream object (defaults to cout)
   *  Note: This function will only print tensors with rank <= 2
   */
  void Print(std::ostream &__os = std::cout) const {
    if (_rank == 0)
      __os << "{Rank = 0 " << _data[0] << " }";
    else if (_rank == 1) {
      __os << "{Rank = 1 ( ";
      for (int mu = 0; mu < 3; mu++)
        __os << _data[mu] << ",";
      __os << _data[3] << ") } ";
    } else if (_rank == 2) {
      int index;
      __os << "{Rank = 2 ";
      for (int mu = 0; mu < 4; mu++) {
        __os << "(";
        for (int nu = 0; nu < 4; nu++) {
          index = 4 * nu + mu;
          __os << _data[index];
          if (nu < 3)
            __os << ",";
        }
        __os << ")";
        if (mu < 3)
          __os << ",";
      }
      __os << "}";
    } else {
      std::cout
          << "<Tensor::Print(ostream&)> Error! Can NOT print a Tensor with "
          << " Rank > 2." << std::endl;
    }
  }

  // Getters:

  /// Returns a constant reference to the @a entry element
  inline const _Tp &operator[](int __entry) const { return _data[__entry]; }

  /// Returns a reference to the @a entry element
  inline _Tp &operator[](int __entry) { return _data[__entry]; }

  /** Returns the \f$ (\mu,\nu,\rho,\sigma,\delta,\pi) \f$ element
   * The arguments all default to zero. Thus, for a rank @a R tensor, only
   *  @a R indicies should be specified. For example, Element(1,2,0) will
   *  access the mu = 1, nu = 2, rho = 0 element of a 3rd rank tensor, etc.
   */
  inline const _Tp &Element(int __mu = 0, int __nu = 0, int __rho = 0,
                            int __sig = 0, int __del = 0, int __pi = 0) const {
    int index = (__pi << 10) + (__del << 8) + (__sig << 6) + (__rho << 4) +
                (__nu << 2) + __mu;
    bool valid = index < this->Size();
    if (!valid)
      std::cout << "Error! Attempt to access non-existant Tensor element."
                << std::endl;
    assert(valid);
    return _data[index];
  }

  /** Returns the \f$ (\mu,\nu,\rho,\sigma,\delta,\pi) \f$ element
   * The arguments all default to zero. Thus, for a rank @a R tensor, only
   *  @a R indicies should be specified. For example, Element(1,2,0) will
   *  access the mu = 1, nu = 2, rho = 0 element of a 3rd rank tensor, etc.
   */
  inline _Tp &Element(int __mu = 0, int __nu = 0, int __rho = 0, int __sig = 0,
                      int __del = 0, int __pi = 0) {
    int index = (__pi << 10) + (__del << 8) + (__sig << 6) + (__rho << 4) +
                (__nu << 2) + __mu;
    bool valid = index < this->Size();
    if (!valid)
      std::cout << "Error! Attempt to access non-existant Tensor element."
                << std::endl;
    assert(valid);
    return _data[index];
  }

  /** Returns the \f$ (\mu,\nu,\rho,\sigma,\delta,\pi) \f$ element
   * See Element for details
   */
  inline const _Tp &operator()(int __mu = 0, int __nu = 0, int __rho = 0,
                               int __sig = 0, int __del = 0,
                               int __pi = 0) const {
    return this->Element(__mu, __nu, __rho, __sig, __del, __pi);
  }

  /** Returns the \f$ (\mu,\nu,\rho,\sigma,\delta,\pi) \f$ element
   * See Element for details
   */
  inline _Tp &operator()(int __mu = 0, int __nu = 0, int __rho = 0,
                         int __sig = 0, int __del = 0, int __pi = 0) {
    return this->Element(__mu, __nu, __rho, __sig, __del, __pi);
  }

  /// Return the element given by @a index (see TensorIndex for details)
  inline _Tp &Element(const TensorIndex &__index) { return _data[__index()]; }

  /// Return the element given by @a index (see TensorIndex for details)
  inline const _Tp &Element(const TensorIndex &__index) const {
    return _data[__index()];
  }

  /// Return the element given by @a index
  inline _Tp &operator()(const TensorIndex &__index) {
    return _data[__index()];
  }

  /// Return the element given by @a index
  inline const _Tp &operator()(const TensorIndex &__index) const {
    return _data[__index()];
  }

  // operators:

  /** Assignment operator
   * Note: Legal if @a Tp = @a T is a legal assignment.
   */
  template <typename T> Tensor<_Tp> &operator=(const Tensor<T> &__tensor) {
    if (!this->RankCheck(__tensor))
      this->SetRank(__tensor.Rank());
    this->Tensor_Base::operator=(__tensor);
    this->_Copy(__tensor);

    return *this;
  }

  /** Assignment operator (rank 0 only)
   * Note: Legal if @a Tp = @a T is a legal assignment.
   */
  template <typename T>
  typename DisableIf<IsTensor(T), Tensor<_Tp> &>::Type operator=(const T &__x) {
    if (_rank != 0)
      std::cout << "Error! Attempt to assign tensor (rank > 0) to a scalar"
                << std::endl;
    assert(_rank == 0);
    _data[0] = __x;

    return *this;
  }

  /// Conversion operator to type @a Tp (valid only for rank 0)
  operator _Tp() const {
    if (_rank != 0) {
      std::cout << "Error! Attempt to convert a tensor (rank != 0) to a scalar."
                << std::endl;
    }
    assert(_rank == 0);
    return _data[0];
  }

  /** Contracts the last index of @a this with the 1st index of @a tensor
   *
   * Returns \f$ R_{\mu_1\mu_2\ldots\nu_1\nu_2\ldots}
   *  = X_{\mu_1\mu_2\ldots\rho} T^{\rho}{}_{\nu_1\nu_2\ldots} \f$
   *
   * Note: Legal if @a Tp * @a T is a legal operation.
   *
   * Return type is Tensor<typename MultType<_Tp,T>::Type> where MultType is
   * the return type of Tp * T (defined in OperationType.h).
   *
   * Example: Define Tensor<float> A(2),B(3), then A*B is a 3rd rank tensor
   * equal to \f$ A_{\mu\nu} B^{\nu}{}_{\pi\delta} \f$
   */
  template <typename T>
  Tensor<typename MultType<_Tp, T>::Type>
  operator*(const Tensor<T> &__tensor) const {
    return this->Contract(__tensor, 1);
  }

  /** Tensor inner product (contracts as many indicies as possible)
   *
   * Returns \f$ R_{\rho\pi\ldots} = X_{\mu_1\mu_2\ldots}
   * T^{\mu_1\mu_2\ldots}{}_{\rho\pi\ldots} \f$ or
   * \f$ X_{\mu_1\mu_2\ldots\rho\pi\ldots} T^{\mu_1\mu_2\ldots} \f$
   * depending on which Tensor has the higher rank.
   *
   * Note: Legal if @a Tp * @a T is a legal operation.
   *
   * Return type is Tensor<typename MultType<_Tp,T>::Type> where MultType is
   * the return type of Tp * T (defined in OperationType.h)
   *
   * Example: Define Tensor<float> A(2),B(3), then (A|B) is a 1st rank tensor
   * equal to \f$ A_{\mu\nu} B^{\mu\nu}{}_{\alpha} \f$
   *
   * Note: If the 2 tensors have the same rank then (A|B) is a rank 0 tensor.
   */
  template <typename T>
  Tensor<typename MultType<_Tp, T>::Type>
  operator|(const Tensor<T> &__tensor) const {
    return this->Contract(__tensor, -1);
  }

  /** Tensor outer product.
   *
   * Returns \f$ R_{\mu_1\mu_2\ldots\nu_1\nu_2\ldots} = X_{\mu_1\mu_2\ldots}
   *  T{\nu_1\nu_2\ldots} \f$
   *
   * Note: Legal if @a Tp * @a T is a legal operation.
   *
   * Return type is Tensor<typename MultType<_Tp,T>::Type> where MultType is
   * the return type of Tp * T (defined in OperationType.h)
   *
   * Example: define Tensor<float> A(2),B(3) then A%B is a 5th rank tensor
   * where A%B = \f$ A_{\mu\nu} B_{\rho\pi\delta} \f$
   */
  template <typename T>
  Tensor<typename MultType<_Tp, T>::Type>
  operator%(const Tensor<T> &__tensor) const {

    int rank = this->_rank + __tensor._rank;
    Tensor<typename MultType<_Tp, T>::Type> ret(rank);

    // if either tensor is rank 0, just return this*__tensor
    if ((_rank == 0) || (__tensor._rank == 0))
      return (*this) * __tensor;
    else { // we actually have to do some work
      TensorIndex index(rank);
      TensorIndex ind1(this->_rank);
      TensorIndex ind2(__tensor._rank);

      int size1 = ind1.Size();
      //    int size2 = ind2.Size();

      while (index.IsValid()) { // loop over ret's elements

        // set up this' indicies
        for (int i = 0; i < size1; i++)
          ind1.SetIndex(i, index[i]);
        // set up _tensor's indicies
        for (int i = size1; i < index.Size(); i++)
          ind2.SetIndex(i - size1, index[i]);

        // set element to product of this and __tensor elements
        ret(index) = (this->Element(ind1)) * (__tensor(ind2));
        index++;
      }
    }
    return ret;
  }

  /** Returns \f$ R_{\mu\nu\ldots} = X_{\mu\nu\ldots} \times x \f$
   * Note: Legal if @a Tp * @a T is a legal operation.
   */
  template <typename T>
  typename EnableIf<IsScalar(T), Tensor<typename MultType<_Tp, T>::Type>>::Type
  operator*(const T &__x) const {
    Tensor<typename MultType<_Tp, T>::Type> ret(_rank);
    int size = this->Size();
    for (int i = 0; i < size; i++)
      ret[i] = _data[i] * __x;
    return ret;
  }

  /** Returns \f$ R_{\mu\nu\ldots} = X_{\mu\nu\ldots} / x \f$
   * Note: Legal if @a Tp / @a T is a legal operation.
   */
  template <typename T>
  typename EnableIf<IsScalar(T), Tensor<typename DivType<_Tp, T>::Type>>::Type
  operator/(const T &__x) const {
    Tensor<typename DivType<_Tp, T>::Type> ret(_rank);
    int size = this->Size();
    for (int i = 0; i < size; i++)
      ret[i] = _data[i] / __x;
    return ret;
  }

  /** Returns \f$ R_{\mu\nu\ldots} = X_{\mu\nu\ldots} + T_{\mu\nu\ldots} \f$
   * Note: Legal if @a Tp + @a T is a legal operation.
   */
  template <typename T>
  Tensor<typename AddType<_Tp, T>::Type>
  operator+(const Tensor<T> &__tensor) const {
    if (this->Rank() != __tensor.Rank())
      std::cout << "Error! Attempt to add tensors w/ different ranks."
                << std::endl;
    assert(_rank == __tensor._rank);
    Tensor<typename AddType<_Tp, T>::Type> ret(_rank);
    int size = this->Size();
    for (int i = 0; i < size; i++)
      ret._data[i] = _data[i] + __tensor._data[i];

    return ret;
  }

  /** Returns \f$ R_{\mu\nu\ldots} = X_{\mu\nu\ldots} - T_{\mu\nu\ldots} \f$
   * Note: Legal if @a Tp - @a T is a legal operation.
   */
  template <typename T>
  Tensor<typename SubType<_Tp, T>::Type>
  operator-(const Tensor<T> &__tensor) const {
    if (this->Rank() != __tensor.Rank())
      std::cout << "Error! Attempt to subtract tensors w/ different ranks."
                << std::endl;
    assert(_rank == __tensor._rank);
    Tensor<typename SubType<_Tp, T>::Type> ret(_rank);
    int size = this->Size();
    for (int i = 0; i < size; i++)
      ret._data[i] = _data[i] - __tensor._data[i];

    return ret;
  }

  /** Rank 0 tensor + scalar
   * Note: Legal if @a Tp + @a T is a legal operation.
   */
  template <typename T>
  typename EnableIf<IsScalar(T), Tensor<typename AddType<_Tp, T>::Type>>::Type
  operator+(const T &__x) const {
    if (_rank != 0)
      std::cout << "Error! Attempt to add tensor (rank > 0) to a scalar"
                << std::endl;
    assert(_rank == 0);
    Tensor<typename AddType<_Tp, T>::Type> ret(0);
    ret() = _data[0] + __x;
    return ret;
  }

  /** Rank 0 tensor - scalar
   * Note: Legal if @a Tp - @a T is a legal operation.
   */
  template <typename T>
  typename EnableIf<IsScalar(T), Tensor<typename SubType<_Tp, T>::Type>>::Type
  operator-(const T &__x) const {
    if (_rank != 0) {
      std::cout << "Error! Attempt to subtract tensor (rank > 0) from a scalar"
                << std::endl;
    }
    assert(_rank == 0);
    Tensor<typename SubType<_Tp, T>::Type> ret(0);
    ret() = _data[0] - __x;
    return ret;
  }

  /// Set @a this = @a this * @a tensor
  template <typename T> Tensor<_Tp> &operator*=(const Tensor<T> &__tensor) {
    (*this) = (*this) * __tensor;
    return *this;
  }

  /// Set @a this = @a this * @a x
  template <typename T>
  typename EnableIf<IsScalar(T), Tensor<_Tp> &>::Type operator*=(const T &__x) {
    (*this) = (*this) * __x;
    return *this;
  }

  /// Set @a this = @a this / @a x
  template <typename T>
  typename EnableIf<IsScalar(T), Tensor<_Tp> &>::Type operator/=(const T &__x) {
    (*this) = (*this) / __x;
    return *this;
  }

  /// Set @a this = @a this + @a tensor
  template <typename T> Tensor<_Tp> &operator+=(const Tensor<T> &__tensor) {
    (*this) = (*this) + __tensor;
    return *this;
  }

  /// Set @a this = @a this - @a tensor
  template <typename T> Tensor<_Tp> &operator-=(const Tensor<T> &__tensor) {
    (*this) = (*this) - __tensor;
    return *this;
  }

  /** This operator shifts the indicies to right @a shift places.
   *
   * Example: define Tensor<float> A(4), then A>>2 returns the tensor
   * \f$ A_{\rho\pi\mu\nu} \f$ if \f$ A = A_{\mu\nu\rho\pi} \f$
   *
   * Note: just returns the tensor if rank < 2
   */
  Tensor operator>>(int __shift) const;

  /** This operator shifts the indicies to the left @a shift places.
   * See operator>> for details.
   */
  Tensor operator<<(int __shift) const;

  /** Comparison operator
   * Requires ranks be the same and each element return @a false under !=
   */
  inline bool operator==(const Tensor<_Tp> &__tensor) const {
    if (!this->RankCheck(__tensor))
      return false;
    int size = this->Size();
    for (int i = 0; i < size; i++)
      if (_data[i] != __tensor._data[i])
        return false;

    return true;
  }

  /// Comparison operator (see operator== for details)
  inline bool operator!=(const Tensor<_Tp> &__tensor) const {
    return !(*this == __tensor);
  }

  // a few extra contraction functions:

  /** Contract 2 tensors.
   *  @param tensor Tensor to contract with @a this
   *  @param num_indicies Number of indicies to contract (defaults to all)
   *
   * Returns \f$ R_{\mu_1\mu_2\ldots\nu_1\nu_2\ldots}
   *  = X_{\mu_1\mu_2\ldots\alpha_1\alpha_1\ldots\alpha_n}
   *  T^{\alpha_1\alpha_2\ldots\alpha_n}{}_{\nu_1\nu_2\ldots} \f$
   *
   * Note: Legal if @a Tp * @a T is a legal operation.
   *
   * Return type is Tensor<typename MultType<_Tp,T>::Type> where MultType is
   * the return type of Tp * T (defined in OperationType.h)
   *
   * Example: define Tensor<float> A(2),B(3), then A.Contract(B) is a 1st rank
   *  tensor equal to \f$ A_{\mu\nu} B^{\mu\nu\delta} \f$, and
   *  A.Contract(B,1) is a 3rd rank tensor equal to
   * \f$ A^{\mu}{}_{\alpha} B^{\alpha\nu\delta} \f$, etc...
   *
   */
  template <typename T>
  Tensor<typename MultType<_Tp, T>::Type>
  Contract(const Tensor<T> &__tensor, int __num_indicies = -1) const {

    int ind1max, ind2st, sumsize, nterm, rank;
    Tensor<typename MultType<_Tp, T>::Type> ret;
    //  MetricTensor g;
    Tensor<double> g(2); // Metric tensor
    g.Element(0, 0) = 1.0;
    g.Element(1, 1) = -1.0;
    g.Element(2, 2) = -1.0;
    g.Element(3, 3) = -1.0;

    double gFactors;
    typename MultType<_Tp, T>::Type element;
    ind1max = 0;
    ind2st = 0;
    sumsize = 0;

    if (_rank ==
        0) { // if this is rank 0, just multiply __tensor by this' value
      ret.SetRank(__tensor._rank);
      ret = _data[0] * __tensor;
      return ret;
    } else if (__tensor._rank == 0) { //__tensor is rank 0...
      ret.SetRank(_rank);
      ret = (*this) * __tensor._data[0];
      return ret;
    }

    if (__num_indicies > _rank || __num_indicies > __tensor._rank) {
      std::cout << "<Tensor::Contract> Error! Can't contract " << __num_indicies
                << " between a " << _rank << " rank and a " << __tensor._rank
                << " rank tensor." << std::endl;
      abort();
    }

    if (__num_indicies > 0)
      rank = _rank + __tensor._rank - 2 * __num_indicies;
    else
      rank = abs(_rank - __tensor._rank);

    // set ret's rank and create the summed TensorIndex
    ret.SetRank(rank);
    TensorIndex indSummed;

    if (__num_indicies < 0) {
      // check to see which tensor has the smaller rank (calculate how many
      // summed indicies are needed)
      if (_rank <= __tensor._rank)
        sumsize = _rank;
      else
        sumsize = __tensor._rank;
    } else
      sumsize = __num_indicies;
    indSummed.Resize(2 * sumsize);

    TensorIndex ind1(this->_rank);
    TensorIndex ind2(__tensor._rank);

    int size1 = ind1.Size();
    int size2 = ind2.Size();

    if (rank > 0) { // rank > 0, we need to do all the loops
      TensorIndex index(rank);

      while (index.IsValid()) { // loop over ret's elements

        // check to see if this will have any free indicies
        if ((size1 - sumsize) > 0)
          ind1max = size1 - sumsize;
        else
          ind1max = 0;

        // set this index (except last ??(number of summed indicies) indicies)
        for (int i = 0; i < ind1max; i++)
          ind1.SetIndex(i, index[i]);

        // set __tensor index (except 1st ??(summed indicies) indicies)
        ind2st = sumsize;
        for (int i = ind2st; i < size2; i++)
          ind2.SetIndex(i, index[ind1max + (i - ind2st)]);

        nterm = 0;
        while (indSummed.IsValid()) { // loop over summed indicies

          gFactors = g(indSummed[0], indSummed[0 + indSummed.Size() / 2]);
          // get the needed amount of metric tensor factors
          for (int i = 1; i < indSummed.Size() / 2; i++) {
            gFactors *= g(indSummed[i], indSummed[i + indSummed.Size() / 2]);
          }
          if (gFactors != 0.0) {
            nterm++;
            // set up last ?? this and 1st ?? __tensor indicies
            for (int i = ind1max; i < size1; i++)
              ind1.SetIndex(i, indSummed[i - ind1max]);
            for (int i = size1 - ind1max; i < indSummed.Size(); i++)
              ind2.SetIndex(i - (size1 - ind1max), indSummed[i]);

            // multiply the metric tensor factor by this and __tensor elements
            element = (this->Element(ind1)) * (__tensor(ind2)) * gFactors;

            // add to each element this*__tensor*g*g...*g with correct # of g's
            if (nterm == 1)
              ret(index) = element;
            else
              ret(index) += element;
          }
          ++indSummed;
        }
        // reset the summed indicies, step up index to next ret element
        indSummed.Reset();
        ++index;
      }
    } else { // both are same rank tensors (R is rank 0)
      nterm = 0;

      // loop over summed indicies (only loop needed in this case)
      while (indSummed.IsValid()) {
        gFactors = g(indSummed[0], indSummed[0 + indSummed.Size() / 2]);

        // get the needed amount of metric tensor factors
        for (int i = 1; i < indSummed.Size() / 2; i++)
          gFactors *= g(indSummed[i], indSummed[i + indSummed.Size() / 2]);

        if (gFactors != 0.0) {
          nterm++;
          for (int i = 0; i < indSummed.Size() / 2; i++) {
            ind1.SetIndex(i, indSummed[i]);
            ind2.SetIndex(i, indSummed[i + indSummed.Size() / 2]);
          }
          element = (this->Element(ind1)) * (__tensor(ind2)) * gFactors;
          if (nterm == 1)
            ret() = element;
          else
            ret() += element;
        }
        indSummed++;
      }
    }
    return ret;
  }

  // some miscelaneous functions:

  /** Permutes the indicies specified by @a mu and @a nu.
   *
   * Example: define Tensor<float> A(3), then A.Permute(0,2) will permute the
   * 1st and 3rd indicies returning \f$ A_{\rho\nu\mu} \f$ if
   * \f$ A = A_{\mu\nu\rho} \f$ (C indexing, starts from zero)
   *
   * Note: just returns the tensor if rank < 2 or either mu or nu >= rank
   */
  Tensor Permute(int __mu, int __nu) const;

  /** Reorder the indicies of the tensor given by @a order.
   *
   * Example: If order = (0,2,1), then Tensor<float> A(3) A.Order(order)
   * returns \f$ A_{\mu\rho\nu} \f$ if \f$ A = A_{\mu\nu\rho} \f$.
   */
  Tensor Order(const TensorIndex &__order) const;

  /** Returns the symmetric tensor built from the current tensor.
   *
   * Example: define Tensor<float> A(3), if \f$ A = A_{\mu\nu\rho} \f$ then
   * A.Symmetric() = \f$ (A_{\mu\nu\rho} + A_{\mu\rho\nu} + A_{\nu\mu\rho}
   *               + A_{\nu\rho\mu} + A_{\rho\mu\nu} + A_{\rho\nu\mu})/6.0 \f$
   */
  Tensor Symmetric() const;

  /** Returns the anti-symmetric Tensor built from the current Tensor.
   *
   * Example: define Tensor<float> A(3), if \f$ A = A_{\mu\nu\rho} \f$ then
   * A.AntiSymmetric() = \f$ (A_{\mu\nu\rho} - A_{\mu\rho\nu} - A_{\nu\mu\rho}
   *               + A_{\nu\rho\mu} + A_{\rho\mu\nu} - A_{\rho\nu\mu})/6.0 \f$
   */
  Tensor AntiSymmetric() const;

  /** Returns the complex conjugate of @a this tensor
   * Note: Legal if conj(_Tp) exists
   *  (see c++-template-utils/TemplateUtilFuncs.h)
   */
  Tensor<_Tp> Conjugate() const {
    Tensor<_Tp> ret(_rank);
    int size = this->Size();
    for (int i = 0; i < size; i++)
      ret[i] = conj(_data[i]);

    return ret;
  }

  /// Return the magnitude squared (\f$ X_{\mu\nu...} X^{\mu\nu...} \f$)
  inline _Tp Mag2() const { return ((*this) | (*this)).Element(); }

  /// Tensor outer product (see operator% for details)
  template <typename T>
  Tensor<typename MultType<_Tp, T>::Type>
  OutterProduct(const Tensor<T> &__tensor) const {
    return (*this) % __tensor;
  }

  /// Tensor inner product (see operator| for details)
  template <typename T>
  Tensor<typename MultType<_Tp, T>::Type>
  InnerProduct(const Tensor<T> &__tensor) const {
    return ((*this) | __tensor);
  }

  /** Lorentz transform the tensor.
   *
   * @param lt A Lorentz transformation tensor (\f$ \Lambda^{\mu}{}_{\nu}\f$)
   *
   *  This function performs a Lorentz transformation on @a this tensor given
   *  by \f$ X_{\mu_1 \mu_2 ...} \Lambda^{\mu_1}{}_{\nu_1}
   *  \Lambda^{\mu_2}{}_{\nu_2} ... \f$
   */
  void Transform(const Tensor<double> &__lt) {

    if (__lt.Rank() != 2)
      std::cout << "Error! Lorentz transformation tensor NOT rank 2."
                << std::endl;
    assert(__lt.Rank() == 2);
    int rank = this->Rank();
    if (rank > 0) { // if rank 0 no transformation needed
      TensorIndex index(rank);
      TensorIndex indSummed(rank);
      int nterm;
      double lamFact;
      // make a copy
      Tensor<_Tp> copy(*this);

      while (index.IsValid()) { // loop over elements of this tensor
        nterm = 0;
        while (indSummed.IsValid()) {
          // get the appropriate number of Lambda_mu_nu factors
          lamFact = __lt(index[0], indSummed[0]);
          for (int i = 1; i < rank; i++)
            lamFact *= __lt(index[i], indSummed[i]);
          if (lamFact != 0.0) {
            nterm++;
            // add to each element this*Lambda*Lambda*...*Lambda
            if (nterm == 1)
              (*this)(index) = lamFact * (copy(indSummed));
            else
              (*this)(index) += lamFact * (copy(indSummed));
          }
          ++indSummed;
        }
        // reset summed indicies, step up index to next element
        indSummed.Reset();
        ++index;
      }
    }
  }
};

//_____________________________________________________________________________
//
// Specifications of functions in TemplateUtilFuncs.h for the Tensor class
//_____________________________________________________________________________

/// Returns a Tensor = T.Zero() (of the same rank as @a tensor)
template <typename T> inline Tensor<T> zero(const Tensor<T> &__tensor) {
  Tensor<T> ret(__tensor.Rank());
  ret.Zero();
  return ret;
}
//_____________________________________________________________________________

/// Same as Tensor::Conjugate
template <typename T> inline Tensor<T> conj(const Tensor<T> &__tensor) {
  return __tensor.Conjugate();
}
//_____________________________________________________________________________

/// Returns a rank 0 tensor with value unity(_Tp)
template <typename T> inline Tensor<T> unity(const Tensor<T> &__tensor) {
  Tensor<T> ret(0);
  T var_type;
  ret() = unity(var_type);
  return ret;
}
//_____________________________________________________________________________

/// Scalar * tensor (see Tensor::operator*)
template <typename T1, typename T2>
typename EnableIf<IsScalar(T1), Tensor<typename MultType<T1, T2>::Type>>::Type
operator*(const T1 &__x, const Tensor<T2> &__tensor) {
  Tensor<typename MultType<T1, T2>::Type> ret(__tensor.Rank());
  int size = __tensor.Size();
  for (int i = 0; i < size; i++)
    ret[i] = __x * __tensor[i];
  return ret;
}
//_____________________________________________________________________________

/// Scalar + rank 0 tensor (see Tensor::operator+)
template <typename T1, typename T2>
typename EnableIf<IsScalar(T1), Tensor<typename AddType<T1, T2>::Type>>::Type
operator+(const T1 &__x, const Tensor<T2> &__tensor) {
  return __tensor + __x;
}
//_____________________________________________________________________________

/// Scalar - rank 0 tensor (see Tensor::operator-)
template <typename T1, typename T2>
typename EnableIf<IsScalar(T1), Tensor<typename SubType<T1, T2>::Type>>::Type
operator-(const T1 &__x, const Tensor<T2> &__tensor) {
  return (__tensor - __x) * -1.;
}
//_____________________________________________________________________________

/// ostream operator for the Tensor class
template <typename T>
inline std::ostream &operator<<(std::ostream &__os, const Tensor<T> &__tensor) {
  __tensor.Print(__os);
  return __os;
}
//_____________________________________________________________________________

/// Return the real part of the tensor
template <typename T> Tensor<T> Real(const Tensor<std::complex<T>> &__tensor) {
  Tensor<T> real(__tensor.Rank());
  for (int i = 0; i < __tensor.Size(); i++)
    real[i] = __tensor[i].real();
  return real;
}
//_____________________________________________________________________________
}
}
#endif /* _Tensor_H */
