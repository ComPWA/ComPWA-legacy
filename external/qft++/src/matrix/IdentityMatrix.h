// IdentityMatrix class definition file -*- C++ -*-
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
#ifndef _IdentityMatrix_H
#define _IdentityMatrix_H

#include "Matrix.h"
#include "Matrix.tcc"
//_____________________________________________________________________________
/** @file IdentityMatrix.h
 *  @brief IdentityMatrix template class definition file.
 */
//_____________________________________________________________________________
/** @class IdentityMatrix
 *  @author Mike Williams 
 *
 *  @brief The identity matrix for Matrix<_Tp>.
 *
 * A square diagnol Matrix. The diagnol elements are set to the value
 * of unity(_Tp). For example, the diagnol elements of IdentiyMatrix<float> are
 * 1.0, while for IdentiyMatrix<complex<float> > they're complex<float>(1.,0.).
 * For more info on unity, see TemplateUtilFuncs.h.
 */
//_____________________________________________________________________________

template <typename _Tp> class IdentityMatrix : public Matrix<_Tp> {

 public:

  /** Constructor
   * @param dim Number of dimensions 
   */
  IdentityMatrix(int __dim):Matrix<_Tp>::Matrix(__dim,__dim){
    _Tp var_type = _Tp();
    for(int i = 0; i < __dim; i++) this->Element(i,i) = unity(var_type);   
  }
    
  /** Destructor */
  ~IdentityMatrix(){}  
};
//_____________________________________________________________________________

#endif /* _IdentityMatrix_H */
