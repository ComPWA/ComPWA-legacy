// Matrix template class source file. -*- C++ -*-
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
#ifndef _Matrix_TCC
#define _Matrix_TCC
//_____________________________________________________________________________
/** @file Matrix.tcc
 *  @brief Matrix template class source file.
 */
//_____________________________________________________________________________

template <typename _Tp> template <typename V> 
Matrix<_Tp>& Matrix<_Tp>::operator=(const Matrix<V> &__matrix) {
  if(!(this->SizeCheck(__matrix))){
    cout << "Error! Attempt to use operator= on matricies w/ different sizes"
	 << endl;
  }
  this->_SizeAssert(__matrix);
  this->Matrix_Base::operator=(__matrix);
  this->_Copy(__matrix);

  return *this;
}
//_____________________________________________________________________________

template <typename _Tp> template <typename V> 
Matrix<typename AddType<_Tp,V>::Type> 
Matrix<_Tp>::operator+(const Matrix<V> &__matrix) const {
  if(!this->SizeCheck(__matrix)){
    cout << "Error! Attempt to add 2 matricies w/ different sizes." << endl;
  }
  this->_SizeAssert(__matrix);
  
  Matrix<typename AddType<_Tp,V>::Type> ret(_num_rows,_num_cols);
  int size = this->Size();
  for(int i = 0; i < size; i++) 
    ret._data[i] = this->_data[i] + __matrix._data[i];
  
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp> template <typename V> 
Matrix<typename SubType<_Tp,V>::Type> 
Matrix<_Tp>::operator-(const Matrix<V> &__matrix) const {
  if(!this->SizeCheck(__matrix)){
    cout << "Error! Attempt to add 2 matricies w/ different sizes." << endl;
  }
  this->_SizeAssert(__matrix);
  
  Matrix<typename SubType<_Tp,V>::Type> ret(_num_rows,_num_cols);
  int size = this->Size();
  for(int i = 0; i < size; i++) 
    ret._data[i] = this->_data[i] - __matrix._data[i];
  
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp> template <typename V> 
Matrix<typename MultType<_Tp,V>::Type> 
Matrix<_Tp>::operator*(const Matrix<V> &__matrix) const {
  if(this->NumCols() != __matrix.NumRows()){
    cout << "Error! Attempt to multiply matricies with incompatible sizes."
	 << endl;
  }
  assert(this->NumCols() == __matrix.NumRows());
  int matrix_cols = __matrix.NumCols();
  Matrix<typename MultType<_Tp,V>::Type> ret(_num_rows,matrix_cols);
  
  for(int i = 0; i < _num_rows; i++){
    for(int j = 0; j < matrix_cols; j++){
      ret(i,j) = this->Element(i,0) * __matrix(0,j);
      for(int k = 1; k < _num_cols; k++) 
	ret(i,j) += this->Element(i,k) * __matrix(k,j);
    }
  }
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp> 
bool Matrix<_Tp>::operator==(const Matrix<_Tp> &__matrix) const {
  if(!this->SizeCheck(__matrix)) return false;
  for(int i = 0; i < _num_rows; i++){
    for(int j = 0; j < _num_rows; j++){
      if(this->Element(i,j) != __matrix(i,j)) return false;
    }
  }
  return true;
}
//_____________________________________________________________________________

template <typename _Tp> _Tp Matrix<_Tp>::Trace() const { 
  if(_num_rows != _num_cols) 
    cout << "Error! Attempt to get trace of non-square matrix." << endl;
  assert(_num_rows == _num_cols);
  _Tp ret;
  ret = this->Element(0,0);
  for(int i = 1; i < _num_rows; i++) ret += this->Element(i,i);

  return ret;
}
//_____________________________________________________________________________

template <typename _Tp>  Matrix<_Tp>  Matrix<_Tp>::Transpose() const {  
  Matrix<_Tp> ret(_num_cols,_num_rows);
  for(int i = 0; i < _num_rows; i++){
    for(int j = 0; j < _num_cols; j++){
      ret(j,i) = this->Element(i,j);
    }
  }
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp> template <typename V> 
Matrix<typename MultType<_Tp,V>::Type> 
Matrix<_Tp>::operator%(const Matrix<V> &__matrix) const {
  if(this->NumCols() != __matrix.NumRows()){
    cout << "Error! Attempt to multiply matricies with incompatible sizes."
	 << endl;
  }
  assert(this->NumCols() == __matrix.NumRows());
  int matrix_cols = __matrix.NumCols();
  Matrix<typename MultType<_Tp,V>::Type> ret(_num_rows,matrix_cols);
  
  for(int i = 0; i < _num_rows; i++){
    for(int j = 0; j < matrix_cols; j++){
      ret(i,j) = this->Element(i,0) % __matrix(0,j);
      for(int k = 1; k < _num_cols; k++) 
	ret(i,j) += this->Element(i,k) % __matrix(k,j);
    }
  }
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp> template <typename V> 
Matrix<typename MultType<_Tp,V>::Type> 
Matrix<_Tp>::operator|(const Matrix<V> &__matrix) const {
  if(this->NumCols() != __matrix.NumRows()){
    cout << "Error! Attempt to multiply matricies with incompatible sizes."
	 << endl;
  }
  assert(this->NumCols() == __matrix.NumRows());
  int matrix_cols = __matrix.NumCols();
  Matrix<typename MultType<_Tp,V>::Type> ret(_num_rows,matrix_cols);
  
  for(int i = 0; i < _num_rows; i++){
    for(int j = 0; j < matrix_cols; j++){
      ret(i,j) = this->Element(i,0) | __matrix(0,j);
      for(int k = 1; k < _num_cols; k++) 
	ret(i,j) += this->Element(i,k) | __matrix(k,j);
    }
  }
  return ret;
}
//_____________________________________________________________________________

#endif /* _Matrix_TCC */
