// Matrix template class definition file. -*- C++ -*-
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
#ifndef _Matrix_H
#define _Matrix_H
//_____________________________________________________________________________
// Standard C++ Headers:
#include <iostream>
#include <cassert>
#include <vector>
#include <complex>
#include <cstdlib>
// Local Headers:
#include "../../include/c++-template-utils.h"
#include "Matrix_Base.h"
//_____________________________________________________________________________
/** @file Matrix.h
 *  @brief Matrix template class definition file.
 */
//_____________________________________________________________________________

using namespace std;
//_____________________________________________________________________________
/** @class Matrix
 *  @author Mike Williams
 *
 *  @brief General template class for handling matricies and matrix operations.
 *
 * Matrix is a template class for handling matricies and matrix operations.
 * A Matrix object can be any size (n x m) and store any type that can be 
 * stored in a C++ container class. This class has been written to be as
 * general and flexible as possible. Parameter passing has been optimized
 * using the Type class. Matrix has been designed so that tow instantiations,
 * Matrix<_A> and Matrix<_B>, are completely compatible as long as _A and _B 
 * are compatible (eg. _A*_B,_A+_B,...are defined).  Return types of 
 * <em>mixed type</em> matrix operations are determined by the OperationType
 * template classes (MultType,etc...).
 *
 * <b>Example Usage </b>
 * 
 * \include Matrix.ex
 */
//_____________________________________________________________________________

template <typename _Tp> class Matrix : public Matrix_Base {

private:
  // attributes:
  vector<_Tp> _data; ///< Matrix elements (type _Tp)

  // private functions:

  /// Copy @a matrix elements to this matrix
  template<typename V> void _Copy(const Matrix<V> &__matrix){
    int size = __matrix.Size();
    _data.resize(size);
    for(int i = 0; i < size; i++) _data[i] = __matrix._data[i];
  }

protected:
  // friends:
  template <typename V> friend class Matrix;

public:

  // create/copy/destroy:

  /** Default Constructor */
  Matrix() : Matrix_Base(){}

  /** Constructor
   * @param nrows Number of rows
   *  @param ncols Number of columns
   */
  Matrix(int __nrows,int __ncols) : Matrix_Base(__nrows,__ncols){
    _data.resize(__nrows*__ncols);
  }

  /** Constructor
   * @param nrows Number of rows
   *  @param ncols Number of columns
   *  @param init  Initial value of Matrix elements
   */
  Matrix(int __nrows,int __ncols,typename Type<_Tp>::ParamType __init)
    : Matrix_Base(__nrows,__ncols){
    _data.resize(__nrows*__ncols,__init);
  }

  /** Constructor
   * @param vv 2-D vector to convert to a Matrix 
   */
  template<typename V> Matrix(const vector<vector<V> > &__vv) : Matrix_Base(){
    _num_rows = (int)__vv.size();
    if(_num_rows == 0) _num_cols = 0; 
    else _num_cols = (int)__vv[0].size();
    _data.resize(_num_rows*_num_cols);

    for(int i = 0; i < _num_rows; i++){
      for(int j = 0; j < _num_cols; j++) (*this)(i,j) = __vv[i][j];
    }
  } 

  /// Copy Constructor
  template<typename V> Matrix(const Matrix<V> &__matrix):Matrix_Base(__matrix){
    this->_Copy(__matrix);
  }

  /** Destructor */
  virtual ~Matrix(){}

  // operators:

  /** Assignment operator
   *
   *  Note: Legal if @a Tp = @a T is a legal assignment. 
   *
   *  Note: Asserts @a matrix and @a this be the same size.
   */
  template <typename V> Matrix<_Tp>& operator=(const Matrix<V> &__matrix);
   
  /// Assignment operator (for 1 x 1 matrix only) 
  Matrix<_Tp>& operator=(typename Type<_Tp>::ParamType __val){
    assert(this->NumRows() == 1 && this->NumCols() == 1);
    _data[0] = __val;
    
    return *this;
  }

  /// Conversion operator to type @a Tp (valid only for 1x1)
  operator _Tp () const {
    if(_num_rows != 1 || _num_cols != 1) {
      cout << "Error! Attempt to convert a matrix (not 1x1) to a scalar." 
	   << endl;
    }
    assert(_num_rows == 1 || _num_cols == 1);
    return _data[0];
  }
  
  /** Matrix addition.
   *
   *  Note: Legal if @a Tp + @a V is a legal operation. 
   *
   *  Note: Asserts @a matrix and @a this be the same size.
   */
  template <typename V> Matrix<typename AddType<_Tp,V>::Type> 
  operator+(const Matrix<V> &__matrix) const;

  /** Matrix subtraction.
   *
   *  Note: Legal if @a Tp - @a V is a legal operation. 
   *
   *  Note: Asserts @a matrix and @a this be the same size.
   */
  template <typename V> Matrix<typename SubType<_Tp,V>::Type> 
  operator-(const Matrix<V> &__matrix) const;

  /** Matrix multiplication
   *
   *  Note: Legal if @a Tp * @a V is a legal operation. 
   *
   *  Note: Asserts the number of rows in @a matrix and the number of columns 
   *        in @a this be the same.
   */
  template <typename V> Matrix<typename MultType<_Tp,V>::Type> 
  operator*(const Matrix<V> &__matrix) const;

  /** Matrix multiplication with % operator
   *
   *  Note: Legal if @a Tp % @a V is a legal operation. 
   *
   *  Note: Asserts the number of rows in @a matrix and the number of columns 
   *        in @a this be the same.
   */
  template <typename V> Matrix<typename MultType<_Tp,V>::Type> 
  operator%(const Matrix<V> &__matrix) const;

  /** Matrix multiplication with | operator
   *
   *  Note: Legal if @a Tp | @a V is a legal operation. 
   *
   *  Note: Asserts the number of rows in @a matrix and the number of columns 
   *        in @a this be the same.
   */
  template <typename V> Matrix<typename MultType<_Tp,V>::Type> 
  operator|(const Matrix<V> &__matrix) const;

  /// \f$ M_{ij} \f$ * x for each matrix element
  template <typename V> 
  typename DisableIf<IsMatrix(V),Matrix<typename MultType<_Tp,V>::Type> >::Type
  operator*(const V &__x) const {
    int size = this->Size();
    Matrix<typename MultType<_Tp,V>::Type> ret(_num_rows,_num_cols);
    for(int i = 0; i < size; i++) ret._data[i] = this->_data[i] * __x;
    return ret;
  }

  /// \f$ M_{ij} \f$ / x for each matrix element
  template <typename V> 
  typename DisableIf<IsMatrix(V),Matrix<typename DivType<_Tp,V>::Type> >::Type
  operator/(const V &__x) const {
    int size = this->Size();
    Matrix<typename DivType<_Tp,V>::Type> ret(_num_rows,_num_cols);
    for(int i = 0; i < size; i++) ret._data[i] = this->_data[i] / __x;
    return ret;
  }
  
  /// \f$ M_{ij} \f$ % x for each matrix element
  template <typename V> 
  typename DisableIf<IsMatrix(V),Matrix<typename MultType<_Tp,V>::Type> >::Type
  operator%(const V &__x) const {
    int size = this->Size();
    Matrix<typename MultType<_Tp,V>::Type> ret(_num_rows,_num_cols);
    for(int i = 0; i < size; i++) ret._data[i] = this->_data[i] % __x;
    return ret;
  }

  /// \f$ M_{ij} \f$ | x for each matrix element
  template <typename V> 
  typename DisableIf<IsMatrix(V),Matrix<typename MultType<_Tp,V>::Type> >::Type
  operator|(const V &__x) const {
    int size = this->Size();
    Matrix<typename MultType<_Tp,V>::Type> ret(_num_rows,_num_cols);
    for(int i = 0; i < size; i++) ret._data[i] = this->_data[i] | __x;
    return ret;
  }

  /// Set @a this = @a this + @a matrix 
  template <typename V> Matrix<_Tp>& operator+=(const Matrix<V> &__matrix){
    (*this) = (*this) + __matrix;
    return *this;
  }

  /// Set @a this = @a this - @a matrix 
  template <typename V> Matrix<_Tp>& operator-=(const Matrix<V> &__matrix){
    (*this) = (*this) - __matrix;
    return *this;
  }

  /// Set @a this = @a this * @a matrix 
  template <typename V> Matrix<_Tp>& operator*=(const Matrix<V> &__matrix){
    (*this) = (*this) * __matrix;
    return *this;
  }

  /// Set @a this = @a this * @a x
  template <typename V> typename DisableIf<IsMatrix(V),Matrix<_Tp>&>::Type 
  operator*=(const V &__x){
    (*this) = (*this) * __x;
    return *this;
  }

  /// Set @a this = @a this / @a x
  template <typename V> typename DisableIf<IsMatrix(V),Matrix<_Tp>&>::Type 
  operator/=(const V &__x){
    (*this) = (*this) / __x;
    return *this;
  }

  /** Comparison operator
   *
   *  @param matrix A Matrix of the same type as @a this
   *
   *  Requires that the matricies be the same size, and that none of their
   *  elements return @a true when passed to the != operator for type @a _Tp.
   */
  bool operator==(const Matrix<_Tp> &__matrix) const;

  /// Comparison operator (see operator== for details)
  inline bool operator!=(const Matrix<_Tp> &__matrix) const {
    return !(*this == __matrix);
  }

  // Getters:

  /// Returns a reference to the (i,j) element of the matrix
  inline _Tp& Element(int __i, int __j){
    if(__i >= _num_rows || __j >= _num_cols)
      cout << "Error! Attempting to access non-existent element." << endl;
    assert(__i < _num_rows && __j < _num_cols);
    return _data[__i*_num_cols + __j];
  }

  /// Returns a constant reference to the (i,j) element of the matrix
  inline const _Tp& Element(int __i, int __j) const {
    if(__i >= _num_rows || __j >= _num_cols)
      cout << "Error! Attempting to access non-existent element." << endl;
    assert(__i < _num_rows && __j < _num_cols);
    return _data[__i*_num_cols + __j];
  }

  /// Returns a reference to the (i,j) element of the matrix
  inline _Tp& operator()(int __i,int __j) {
    return this->Element(__i,__j);
  }

  /// Returns a constant reference to the (i,j) element of the matrix
  inline const _Tp& operator()(int __i,int __j) const {
    return this->Element(__i,__j);
  }

  /// Returns a constant reference to the i'th array element
  inline const _Tp& operator[](int __i) const {
    return _data[__i];
  }

  /// Returns the matrix \f$ M_{jk}(i) \f$ from this matrix
  Matrix<_Tp> operator()(int __i) const {
    Matrix<_Tp> ret(_num_rows,_num_cols);
    for(int j = 0; j < this->Size(); j++) 
      ret._data[j] = this->_data[j].operator()(__i);
    
    return ret;
  }

  
  // Functions:

  /** Returns the trace of the matrix */
  _Tp Trace() const;

  /// Returns the trace of the matrix
  inline _Tp Tr() const {
    return this->Trace();
  }

  /** Returns the transpose of @a this matrix */
  Matrix<_Tp> Transpose() const;

  /// Return the transpose of @a this matrix
  inline Matrix<_Tp> T() const {
    return this->Transpose();
  }

  /** Returns the complex conjugate of @a this matrix
   * Note: Legal if conj(_Tp) exists 
   *  (see c++-template-utils/TemplateUtilFuncs.h) 
   */
  Matrix<_Tp> Conjugate() const {
    Matrix<_Tp> ret(_num_rows,_num_cols);
    for(int i = 0; i < this->Size(); i++) ret._data[i] = conj(this->_data[i]);
    return ret;
  }

  /// Returns the complex conjugate of @a this matrix (see Conjugate()).
  Matrix<_Tp> Conj() const {
    return this->Conjugate();
  }

  /// Returns the adjoint of @a this matrix (see Conjugate() and Transpose())
  Matrix<_Tp> Adjoint() const {
    Matrix<_Tp> ret(_num_rows,_num_cols);
    ret = (this->Transpose()).Conjugate();
    return ret;
  }

  /// Returns the adjoint of @a this matrix (see Conjugate() and Transpose())
  Matrix<_Tp> Dagger() const {
    return this->Adjoint();
  }

  /// Set each element of @a this matrix to zero
  /** Note: Legal if zero(_Tp) is legal 
   *  (see c++-template-utils/TemplateUtilFuncs.h)
   */
  void Zero() {
    for(int i = 0; i < this->Size(); i++) _data[i] = zero(_data[i]);
  }

  /// Removes all elements from the matrix
  void Clear() {
    if(!_data.empty()) _data.clear();
    _num_rows = 0;
    _num_cols = 0;
  }

  /// Resize the matrix to be @a nrows X @a ncols.
  void Resize(int __nrows,int __ncols) {
    _num_rows = __nrows;
    _num_cols = __ncols;
    _data.resize(__nrows*__ncols);
  }

  /// Print the matrix to @a os
  void Print(std::ostream &__os = std::cout) const {
    for(int i = 0; i < _num_rows; i++){
      for(int j = 0; j < _num_cols; j++) __os << this->Element(i,j) << " ";
      __os << endl;
    }
  }

};
//_____________________________________________________________________________

/// @a x * \f$ M_{ij} \f$ for each matrix element
template<typename T1,typename T2> 
typename DisableIf<IsMatrix(T1),Matrix<typename MultType<T1,T2>::Type> >::Type
operator*(const T1 &__x,const Matrix<T2> &__matrix){
  int rows = __matrix.NumRows();
  int cols = __matrix.NumCols();
  Matrix<typename MultType<T1,T2>::Type> ret(rows,cols);

  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++){
      ret(i,j) = __x * __matrix(i,j);
    }
  }
  return ret;
}  
//_____________________________________________________________________________

/// @a x % \f$ M_{ij} \f$ for each matrix element
template<typename T1,typename T2> 
typename DisableIf<IsMatrix(T1),Matrix<typename MultType<T1,T2>::Type> >::Type
operator%(const T1 &__x,const Matrix<T2> &__matrix){
  int rows = __matrix.NumRows();
  int cols = __matrix.NumCols();
  Matrix<typename MultType<T1,T2>::Type> ret(rows,cols);

  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++){
      ret(i,j) = __x % __matrix(i,j);
    }
  }
  return ret;
}  
//_____________________________________________________________________________

/// @a x | \f$ M_{ij} \f$ for each matrix element
template<typename T1,typename T2> 
typename DisableIf<IsMatrix(T1),Matrix<typename MultType<T1,T2>::Type> >::Type
operator|(const T1 &__x,const Matrix<T2> &__matrix){
  int rows = __matrix.NumRows();
  int cols = __matrix.NumCols();
  Matrix<typename MultType<T1,T2>::Type> ret(rows,cols);

  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++){
      ret(i,j) = __x | __matrix(i,j);
    }
  }
  return ret;
}  
//_____________________________________________________________________________

/// Overload of operator << for the Matrix class
template<typename T> 
ostream& operator<<(ostream &__os,const Matrix<T> &__matrix){
  __matrix.Print(__os);
  return __os;
}
//_____________________________________________________________________________

#include "Matrix.tcc" /* Matrix template class source file */

#endif /* _Matrix_H */
