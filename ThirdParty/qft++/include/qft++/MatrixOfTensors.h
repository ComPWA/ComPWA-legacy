// File for Matrix/Tensor class usage. -*- C++ -*-
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
#ifndef _MatrixOfTensors_H
#define _MatrixOfTensors_H
//_____________________________________________________________________________
#include "../../include/matrix.h"
#include "../../include/tensor.h"

namespace ComPWA {
    namespace QFT {
      
//_____________________________________________________________________________
/** @file MatrixOfTensors.h
 *  @brief Defines functions for Matrix/Tensor class interactions.
 */
//_____________________________________________________________________________
//
// Functions for using Matrix<Tensor<T> >
//_____________________________________________________________________________

/// Set the rank of the Tensor elements of the Matrix
template<typename T> void SetRank(Matrix<Tensor<T> > &__matrix,int __rank){
  for(int i = 0; i < __matrix.NumRows(); i++){
    for(int j = 0; j < __matrix.NumCols(); j++) __matrix(i,j).SetRank(__rank);
  }
}
//_____________________________________________________________________________

/** Returns the \f$ M^{ij}_{\mu\nu\rho\ldots} \f$
 *
 * From a Matrix of Tensors, it returns the Matrix of just the element of 
 * each Tensor returned by T(index).
 */
template<typename T> 
Matrix<T> TensorElement(const Matrix<Tensor<T> > &__matrix,
			const TensorIndex &__index){
  int rows = __matrix.NumRows();
  int cols = __matrix.NumCols();
  Matrix<T> ret(rows,cols);
  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++) ret(i,j) = __matrix(i,j).Element(__index);  
  }
  return ret;
}
//_____________________________________________________________________________

/// Sets \f$ M^{ij}_{\mu\nu\rho\ldots} = mval^{ij} \f$ 
template <typename T>
void SetTensorElement(Matrix<Tensor<T> > &__matrix,
		      const TensorIndex &__index,
		      const Matrix<Tensor<T> > &__mval){
  int rows = __matrix.NumRows();
  int cols = __matrix.NumCols();
  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++) 
      __matrix(i,j).Element(__index) = __mval(i,j).Element();    
  }
}
//_____________________________________________________________________________

/// Boost a Matrix of Tensors by boost vector (bx,by,bz)
template <typename T>
Matrix<Tensor<T> > BoostMatrix(const Matrix<Tensor<T> > &__matrix,
			       double __bx,double __by,double __bz){
  int rows = __matrix.NumRows();
  int cols = __matrix.NumCols();  
  Matrix<Tensor<T> > ret(__matrix);

  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++) ret(i,j).Boost(__bx,__by,__bz);    
  }
  return ret;
}
//_____________________________________________________________________________
    }
}

#endif /* _MatrixOfTensors_H */
