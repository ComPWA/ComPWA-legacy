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
// Local Headers:
#include "../../include/c++-template-utils.h"
#include "Tensor_Base.h"
#include "TensorIndex.h"
//_____________________________________________________________________________
/** @file Tensor.h
 *  @brief Tensor template class definition file.
 */
//_____________________________________________________________________________

using namespace std;
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
  vector<_Tp> _data; ///< Tensor elements (type _Tp)

private:

  // private functions:

  /// Copy @a tensor elements to @a this tensor
  template<typename T> void _Copy(const Tensor<T> &__tensor){
    int size = __tensor.Size();
    _data.resize(size);
    for(int i = 0; i < size; i++) _data[i] = __tensor._data[i];
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
  Tensor(int __rank) : Tensor_Base(__rank),_data(1 << (__rank << 1)){}

  /** Constructor
   * @param rank Rank of the Tensor
   *  @param init  Initial value of Tensor elements
   */
  Tensor(int __rank,typename Type<_Tp>::ParamType __init): 
    Tensor_Base(__rank),_data(1 << (__rank << 1),__init){}

  /// Copy Constructor
  template<typename T> Tensor(const Tensor<T> &__tensor):Tensor_Base(__tensor){
    this->_Copy(__tensor);
  }

  /** Destructor */
  virtual ~Tensor(){}

  // basic functions:
  
  /// Returns the number of elements in the tensor
  inline int Size() const {
    return(1 << (_rank << 1));
  }

  /** Set each element of @a this tensor to zero
   * Note: Legal if zero(_Tp) is legal 
   *  (see c++-template-utils/TemplateUtilFuncs.h)
   */
  inline void Zero() {
    _Tp var_type;
    for(int i = 0; i < this->Size(); i++) _data[i] = zero(var_type);
  }

  /// Removes all elements from the tensor
  void Clear() {
    if(!_data.empty()) _data.clear();
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
  void Boost(double __bx,double __by,double __bz);

  /// Boost the Tensor to the rest frame of the 4-momentum @a p4.
  void Boost(const Tensor<double> &__p4) {
    if(__p4.Rank() != 1) cout << "Error! 4-momentum NOT rank 1." << endl;
    assert(__p4.Rank() == 1);

    this->Boost(-(__p4(1)/__p4(0)),-(__p4(2)/__p4(0)),-(__p4(3)/__p4(0)));
  }

  /// Rotate the tensor using Euler angles \f$ \alpha,\beta,\gamma \f$.
  void Rotate(double __alpha,double __beta,double __gamma);

  /// Rotate about the x-axis
  void RotateX(double __alpha);

  /// Rotate about the y-axis
  void RotateY(double __alpha);

  /// Rotate about the z-axis
  void RotateZ(double __alpha);
  
  /** Send the values of the tensor elements to @a os.
   *
   *  @param os ostream object (defaults to cout)
   *  Note: This function will only print tensors with rank <= 2
   */
  void Print(std::ostream& __os = std::cout) const;

  // Getters:

  /// Returns a constant reference to the @a entry element
  inline const _Tp& operator[](int __entry) const {
    return _data[__entry];
  }

  /// Returns a reference to the @a entry element
  inline _Tp& operator[](int __entry) {
    return _data[__entry];
  }

  /** Returns the \f$ (\mu,\nu,\rho,\sigma,\delta,\pi) \f$ element
   * The arguments all default to zero. Thus, for a rank @a R tensor, only
   *  @a R indicies should be specified. For example, Element(1,2,0) will
   *  access the mu = 1, nu = 2, rho = 0 element of a 3rd rank tensor, etc.
   */
  inline const _Tp& Element(int __mu = 0,int __nu = 0,int __rho = 0,
			    int __sig = 0,int __del = 0,int __pi = 0) const {
    int index = (__pi << 10) + (__del << 8) + (__sig << 6) + (__rho << 4) 
      + (__nu << 2) + __mu;
    bool valid = index < this->Size();
    if(!valid) 
      cout << "Error! Attempt to access non-existant Tensor element." << endl;
    assert(valid);
    return _data[index];
  }

  /** Returns the \f$ (\mu,\nu,\rho,\sigma,\delta,\pi) \f$ element
   * The arguments all default to zero. Thus, for a rank @a R tensor, only
   *  @a R indicies should be specified. For example, Element(1,2,0) will
   *  access the mu = 1, nu = 2, rho = 0 element of a 3rd rank tensor, etc.
   */
  inline _Tp& Element(int __mu = 0,int __nu = 0,int __rho = 0,int __sig = 0,
		      int __del = 0,int __pi = 0) {
    int index = (__pi << 10) + (__del << 8) + (__sig << 6) + (__rho << 4) 
      + (__nu << 2) + __mu;
    bool valid = index < this->Size();
    if(!valid) 
      cout << "Error! Attempt to access non-existant Tensor element." << endl;
    assert(valid);
    return _data[index];
  }

  /** Returns the \f$ (\mu,\nu,\rho,\sigma,\delta,\pi) \f$ element
   * See Element for details 
   */
  inline const _Tp& operator()(int __mu = 0,int __nu = 0,int __rho = 0,
			       int __sig = 0,int __del = 0,int __pi = 0)const{
    return this->Element(__mu,__nu,__rho,__sig,__del,__pi);
  }

  /** Returns the \f$ (\mu,\nu,\rho,\sigma,\delta,\pi) \f$ element
   * See Element for details 
   */
  inline _Tp& operator()(int __mu = 0,int __nu = 0,int __rho = 0,int __sig = 0,
			 int __del = 0,int __pi = 0) {
    return this->Element(__mu,__nu,__rho,__sig,__del,__pi);
  }

  /// Return the element given by @a index (see TensorIndex for details)   
  inline _Tp& Element(const TensorIndex &__index) {
    return _data[__index()];
  }

  /// Return the element given by @a index (see TensorIndex for details)   
  inline const _Tp& Element(const TensorIndex &__index) const {
    return _data[__index()];
  }

  /// Return the element given by @a index
  inline _Tp& operator()(const TensorIndex &__index) {
    return _data[__index()];
  }

  /// Return the element given by @a index
  inline const _Tp& operator()(const TensorIndex &__index) const {
    return _data[__index()];
  }

  // operators:

  /** Assignment operator
   * Note: Legal if @a Tp = @a T is a legal assignment.
   */
  template<typename T> Tensor<_Tp>& operator=(const Tensor<T> &__tensor){
    if(!this->RankCheck(__tensor)) this->SetRank(__tensor.Rank());
    this->Tensor_Base::operator=(__tensor);
    this->_Copy(__tensor);

    return *this;
  }

  /** Assignment operator (rank 0 only)
   * Note: Legal if @a Tp = @a T is a legal assignment.
   */
  template<typename T> 
  typename DisableIf<IsTensor(T),Tensor<_Tp>&>::Type operator=(const T &__x){
    if(_rank != 0) 
      cout << "Error! Attempt to assign tensor (rank > 0) to a scalar" << endl;
    assert(_rank == 0);
    _data[0] = __x;

    return *this;
  }
  
  /// Conversion operator to type @a Tp (valid only for rank 0)
  operator _Tp () const {
    if(_rank != 0) {
      cout << "Error! Attempt to convert a tensor (rank != 0) to a scalar." 
	   << endl;
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
  template<typename T> Tensor<typename MultType<_Tp,T>::Type>
  operator*(const Tensor<T> &__tensor) const {
    return this->Contract(__tensor,1);
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
  template<typename T> Tensor<typename MultType<_Tp,T>::Type>
  operator|(const Tensor<T> &__tensor) const {
    return this->Contract(__tensor,-1);
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
  template<typename T> Tensor<typename MultType<_Tp,T>::Type>
  operator%(const Tensor<T> &__tensor) const;

  /** Returns \f$ R_{\mu\nu\ldots} = X_{\mu\nu\ldots} \times x \f$
   * Note: Legal if @a Tp * @a T is a legal operation.
   */
  template<typename T> 
  typename EnableIf<IsScalar(T),Tensor<typename MultType<_Tp,T>::Type> >::Type 
  operator*(const T &__x) const {
    Tensor<typename MultType<_Tp,T>::Type> ret(_rank);
    int size = this->Size();
    for(int i = 0; i < size; i++) ret[i] = _data[i] * __x;
    return ret;  
  }

  /** Returns \f$ R_{\mu\nu\ldots} = X_{\mu\nu\ldots} / x \f$
   * Note: Legal if @a Tp / @a T is a legal operation.
   */
  template<typename T> 
  typename EnableIf<IsScalar(T),Tensor<typename DivType<_Tp,T>::Type> >::Type 
  operator/(const T &__x) const {
    Tensor<typename DivType<_Tp,T>::Type> ret(_rank);
    int size = this->Size();
    for(int i = 0; i < size; i++) ret[i] = _data[i] / __x;
    return ret;  
  }

  /** Returns \f$ R_{\mu\nu\ldots} = X_{\mu\nu\ldots} + T_{\mu\nu\ldots} \f$
   * Note: Legal if @a Tp + @a T is a legal operation.
   */
  template<typename T> Tensor<typename AddType<_Tp,T>::Type>
  operator+(const Tensor<T> &__tensor) const {
    if(this->Rank() != __tensor.Rank())
      cout << "Error! Attempt to add tensors w/ different ranks." << endl;
    assert(_rank == __tensor._rank);
    Tensor<typename AddType<_Tp,T>::Type> ret(_rank); 
    int size = this->Size();
    for(int i = 0; i < size; i++) ret._data[i] = _data[i] + __tensor._data[i];

    return ret;
  }

  /** Returns \f$ R_{\mu\nu\ldots} = X_{\mu\nu\ldots} - T_{\mu\nu\ldots} \f$
   * Note: Legal if @a Tp - @a T is a legal operation.
   */
  template<typename T> Tensor<typename SubType<_Tp,T>::Type>
  operator-(const Tensor<T> &__tensor) const {
    if(this->Rank() != __tensor.Rank())
      cout << "Error! Attempt to subtract tensors w/ different ranks." << endl;
    assert(_rank == __tensor._rank);
    Tensor<typename SubType<_Tp,T>::Type> ret(_rank); 
    int size = this->Size();
    for(int i = 0; i < size; i++) ret._data[i] = _data[i] - __tensor._data[i];

    return ret;
  }

  /** Rank 0 tensor + scalar
   * Note: Legal if @a Tp + @a T is a legal operation.
   */
  template<typename T> 
  typename EnableIf<IsScalar(T),Tensor<typename AddType<_Tp,T>::Type> >::Type
  operator+(const T &__x) const {
    if(_rank != 0) 
      cout << "Error! Attempt to add tensor (rank > 0) to a scalar" << endl;
    assert(_rank == 0);
    Tensor<typename AddType<_Tp,T>::Type> ret(0);
    ret() = _data[0] + __x;
    return ret;  
  }

  /** Rank 0 tensor - scalar
   * Note: Legal if @a Tp - @a T is a legal operation.
   */
  template<typename T> 
  typename EnableIf<IsScalar(T),Tensor<typename SubType<_Tp,T>::Type> >::Type
  operator-(const T &__x) const {
    if(_rank != 0){ 
      cout << "Error! Attempt to subtract tensor (rank > 0) from a scalar" 
	   << endl;
    }
    assert(_rank == 0);
    Tensor<typename SubType<_Tp,T>::Type> ret(0);
    ret() = _data[0] - __x;
    return ret;  
  }

  /// Set @a this = @a this * @a tensor
  template<typename T> Tensor<_Tp>& operator*=(const Tensor<T> &__tensor){
    (*this) = (*this) * __tensor;
    return *this;
  }

  /// Set @a this = @a this * @a x
  template<typename T> typename EnableIf<IsScalar(T),Tensor<_Tp>&>::Type
  operator*=(const T &__x){
    (*this) = (*this) * __x;
    return *this;
  }

  /// Set @a this = @a this / @a x
  template<typename T> typename EnableIf<IsScalar(T),Tensor<_Tp>&>::Type
  operator/=(const T &__x){
    (*this) = (*this) / __x;
    return *this;
  }

  /// Set @a this = @a this + @a tensor
  template<typename T> Tensor<_Tp>& operator+=(const Tensor<T> &__tensor) {
    (*this) = (*this) + __tensor;
    return *this;
  }
  
  /// Set @a this = @a this - @a tensor
  template<typename T> Tensor<_Tp>& operator-=(const Tensor<T> &__tensor) {
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
    if(!this->RankCheck(__tensor)) return false;
    int size = this->Size();
    for(int i = 0; i < size; i++) 
      if(_data[i] != __tensor._data[i]) return false;

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
  template<typename T> Tensor<typename MultType<_Tp,T>::Type>
  Contract(const Tensor<T> &__tensor,int __num_indicies = -1) const; 

  // some miscelaneous functions:

  /** Permutes the indicies specified by @a mu and @a nu.
   *
   * Example: define Tensor<float> A(3), then A.Permute(0,2) will permute the 
   * 1st and 3rd indicies returning \f$ A_{\rho\nu\mu} \f$ if 
   * \f$ A = A_{\mu\nu\rho} \f$ (C indexing, starts from zero)
   *
   * Note: just returns the tensor if rank < 2 or either mu or nu >= rank
   */
  Tensor Permute(int __mu,int __nu) const;
    
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
    for(int i = 0; i < size; i++) ret[i] = conj(_data[i]);
  
    return ret;
  }

  /// Return the magnitude squared (\f$ X_{\mu\nu...} X^{\mu\nu...} \f$)
  inline _Tp Mag2() const {
    return ((*this)|(*this)).Element();
  }

  /// Tensor outer product (see operator% for details)
  template<typename T> Tensor<typename MultType<_Tp,T>::Type>
  OutterProduct(const Tensor<T> &__tensor) const {
    return (*this) % __tensor;
  }

  /// Tensor inner product (see operator| for details)
  template<typename T> Tensor<typename MultType<_Tp,T>::Type>
  InnerProduct(const Tensor<T> &__tensor) const {
    return ( (*this) | __tensor);
  }

  /** Lorentz transform the tensor.
   *
   * @param lt A Lorentz transformation tensor (\f$ \Lambda^{\mu}{}_{\nu}\f$)
   *
   *  This function performs a Lorentz transformation on @a this tensor given
   *  by \f$ X_{\mu_1 \mu_2 ...} \Lambda^{\mu_1}{}_{\nu_1} 
   *  \Lambda^{\mu_2}{}_{\nu_2} ... \f$
   */
  void Transform(const Tensor<double> &__lt);

};
//_____________________________________________________________________________

/// Scalar * tensor (see Tensor::operator*)
template<typename T1,typename T2> 
typename EnableIf<IsScalar(T1),Tensor<typename MultType<T1,T2>::Type> >::Type
operator*(const T1 &__x,const Tensor<T2> &__tensor){
  Tensor<typename MultType<T1,T2>::Type> ret(__tensor.Rank());
  int size = __tensor.Size();
  for(int i = 0; i < size; i++) ret[i] = __x * __tensor[i];
  return ret;  
}
//_____________________________________________________________________________

/// Scalar + rank 0 tensor (see Tensor::operator+)
template<typename T1,typename T2> 
typename EnableIf<IsScalar(T1),Tensor<typename AddType<T1,T2>::Type> >::Type
operator+(const T1 &__x,const Tensor<T2> &__tensor){
  return __tensor + __x;
}
//_____________________________________________________________________________

/// Scalar - rank 0 tensor (see Tensor::operator-)
template<typename T1,typename T2> 
typename EnableIf<IsScalar(T1),Tensor<typename SubType<T1,T2>::Type> >::Type
operator-(const T1 &__x,const Tensor<T2> &__tensor){
  return (__tensor - __x) * -1.;
}
//_____________________________________________________________________________

/// ostream operator for the Tensor class
template <typename T>
inline std::ostream& operator<<(std::ostream& __os,const Tensor<T> &__tensor){
  __tensor.Print(__os);
  return __os;
}
//_____________________________________________________________________________

/// Return the real part of the tensor
template <typename T>
Tensor<T> Real(const Tensor<complex<T> > &__tensor){
  Tensor<T> real(__tensor.Rank());
  for(int i = 0; i < __tensor.Size(); i++) real[i] = __tensor[i].real();
  return real;
}
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

#endif /* _Tensor_H */
