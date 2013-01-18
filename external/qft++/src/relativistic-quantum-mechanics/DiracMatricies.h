// Definition file for Dirac matricies classes -*- C++ -*-
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
#ifndef _DiracMatricies_H
#define _DiracMatricies_H
//_____________________________________________________________________________
/** @file DiracMatricies.h
 *  @brief Definition file of DiracGamma,DiracSigma,DiracGamma5 and PauliSigma
 */
//_____________________________________________________________________________

#include "../../include/matrix.h"
#include "../../include/tensor.h"

using namespace std;
//_____________________________________________________________________________
/** @class DiracGamma
 *  @author Mike Williams
 *
 *  @brief \f$ \gamma^{\mu} \f$ : The Dirac Gamma matricies.
 *
 * In general, these matrices satisfy \f$ \{\gamma^{\mu},\gamma^{\nu}
 * \} = 2 g^{\mu\nu} \f$. We have chosen to use the following 
 * representation (written here in 2x2 block diagnol form):
 *
 * \f$ \gamma^0 = \left(\begin{array}{cc} I & 0 \\ 0 & -I \end{array}\right)\f$
 *
 * \f$ \gamma^i = \left(\begin{array}{cc} 0 & \sigma^i 
 * \\ -\sigma^i & 0 \end{array}\right)\f$
 *
 * Where \f$ \sigma^i \f$ are the Pauli sigma matricies (see PauliSigma).
 */
//_____________________________________________________________________________

class DiracGamma : public Matrix<Tensor<complex<double> > > {

public:

  /// Constructor
  DiracGamma():Matrix<Tensor<complex<double> > >::Matrix(4,4){
    for(int i = 0; i < 4; i++){
      for(int j = 0; j < 4; j++) this->Element(i,j).SetRank(1);      
    }
    this->Element(0,0).Element(0) = 1.0;
    this->Element(0,2).Element(3) = 1.0;
    this->Element(0,3).Element(1) = 1.0;
    this->Element(0,3).Element(2) = complex<double>(0.,-1.);
    this->Element(1,1).Element(0) = 1.0;
    this->Element(1,2).Element(1) = 1.0;
    this->Element(1,2).Element(2) = complex<double>(0.,1.);
    this->Element(1,3).Element(3) = -1.0;
    this->Element(2,0).Element(3) = -1.0;
    this->Element(2,1).Element(1) = -1.0;
    this->Element(2,1).Element(2) = complex<double>(0.,1.);
    this->Element(2,2).Element(0) = -1.0;
    this->Element(3,0).Element(1) = -1.0;
    this->Element(3,0).Element(2) = complex<double>(0.,-1.);
    this->Element(3,1).Element(3) = 1.0;
    this->Element(3,3).Element(0) = -1.0;
  };

  /// Destructor
  ~DiracGamma(){}

  /** Returns \f$ \gamma^{\mu} \f$
   * Returns the 4 x 4 Matrix constructed from the @a mu element of each of 
   * the Tensor elements stored.
   */
  Matrix<complex<double> > operator()(int __mu) const {
    Matrix<complex<double> > ret(4,4);
    for(int i = 0; i < 4; i++){
      for(int j = 0; j < 4; j++) ret(i,j) = this->Element(i,j).Element(__mu);
    }
    return ret;
  }

  /// Same as operator()
  inline Matrix<complex<double> > operator[](int __mu)const {
    return (*this)(__mu);
  }

};
//_____________________________________________________________________________
/** @class DiracGamma5
 *  @author Mike Williams
 *
 *  @brief \f$ \gamma^5 \f$: The Dirac Gamma matrix.
 *
 * \f$ \gamma^5 = i \gamma^0\gamma^1\gamma^2\gamma^3 \f$ or, written in 
 * 2 x 2 block diagnol form:
 *
 * \f$ \gamma^5 = \left(\begin{array}{cc} 0 & I \\ I & 0 \end{array} \right)\f$
 */
//_____________________________________________________________________________

class DiracGamma5 : public Matrix<double> {

public:

  /// Constructor
  DiracGamma5():Matrix<double>::Matrix(4,4){
    this->Element(0,2) = 1.0;
    this->Element(1,3) = 1.0;
    this->Element(2,0) = 1.0;
    this->Element(3,1) = 1.0;
  }

  /// Destructor
  ~DiracGamma5(){}

};
//_____________________________________________________________________________
/** @class DiracSigma
 *  @brief \f$ \sigma^{\mu\nu} \f$ : The Dirac sigma matricies
 *
 *  \f$ \sigma^{\mu\nu} = \frac{i}{2}[\gamma^{\mu},\gamma^{\nu}] \f$
 */
//_____________________________________________________________________________

class DiracSigma : public Matrix<Tensor<complex<double> > > {

public:

  /// Constructor
  DiracSigma() : Matrix<Tensor<complex<double> > >::Matrix(4,4){
    for(int i = 0; i < this->NumRows(); i++){
      for(int j = 0; j < this->NumCols(); j++){
	this->Element(i,j).SetRank(2);
      }
    }
    IdentityMatrix<double> I(4);
    DiracGamma gamma;
    MetricTensor g;
    *this = complex<double>(0.,1.)*(gamma%gamma - g*I);
  }

  /// Destructor
  ~DiracSigma(){}

  /// Assignment operator
  DiracSigma& operator=(const Matrix<Tensor<complex<double> > > &__mt){
    for(int i = 0; i < this->NumRows(); i++){
      for(int j = 0; j < this->NumCols(); j++){
	this->Element(i,j) = __mt(i,j);
      }
    }
    return *this;
  }

};
//_____________________________________________________________________________

/** @class PauliSigma
 *  @author Mike Williams
 *
 *  @brief \f$ \sigma \f$ : The Pauli sigma matricies.
 *
 * In general, these matricies satisfy\f$\{\sigma^i,\sigma^j\}=2\delta^{ij}\f$.
 * We have chosen the following representation:
 *
 * \f$ \sigma^1 = \left(\begin{array}{cc}0&1\\1&0 \end{array}\right) \f$
 *
 * \f$ \sigma^2 = \left(\begin{array}{cc}0&-i\\i&0 \end{array}\right) \f$
 *
 * \f$ \sigma^3 = \left(\begin{array}{cc}1&0\\0&-1 \end{array}\right) \f$
 *
 * Note: If we define, PauliSigma sig, then sig(0) is the 2x2 identity matrix.
 */
//_____________________________________________________________________________

class PauliSigma : public Matrix<Tensor<complex<double> > > {

public:

  /// Constructor
  PauliSigma():Matrix<Tensor<complex<double> > >::Matrix(2,2){
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++) this->Element(i,j).SetRank(1);      
    }
    this->Element(0,0).Element(0) = 1.;
    this->Element(0,0).Element(3) = 1.;
    this->Element(0,1).Element(1) = 1.;
    this->Element(0,1).Element(2) = complex<double>(0.,-1.);
    this->Element(1,0).Element(1) = 1.;
    this->Element(1,0).Element(2) = complex<double>(0.,1.);
    this->Element(1,1).Element(0) = 1.;
    this->Element(1,1).Element(3) = -1.;
  }

  /// Destructor
  ~PauliSigma(){}

  /// Returns \f$ \sigma^i \f$.
  Matrix<complex<double> > operator()(int __mu) const {
    Matrix<complex<double> > ret(2,2);
    ret(0,0) = this->Element(0,0).Element(__mu);
    ret(1,0) = this->Element(1,0).Element(__mu);
    ret(0,1) = this->Element(0,1).Element(__mu);
    ret(1,1) = this->Element(1,1).Element(__mu);
    return ret;
  }

  /// Same as operator()
  inline Matrix<complex<double> > operator[](int __mu) const {
    return (*this)(__mu);
  }
};
//_____________________________________________________________________________

/** \f$ \bar{matrix} \f$
 *
 * Bar is defined as follows:
 *
 * Spinors (4x1): \f$ \bar{M} = M^{\dagger} \gamma^0 \f$
 *
 * Projectors (4x4): \f$ \bar{M} = \gamma^0 M^{\dagger} \gamma^0 \f$
 *
 * All other matricies: \f$ \bar{M} = M^{\dagger} \f$
 *
 *<!--
 * Define: Matrix<V> spinor(4,1),projector(4,4), then
 *   Bar(spinor) = (spinor.Adjoint())*gamma(0)
 *   Bar(projector) = gamma(0)*(projector.Adjoint())*gamma(0)
 *
 * For other Matricies it just returns the Adjoint.
 * -->
 */
template <typename T>
Matrix<typename MultType<T,complex<double> >::Type> 
Bar(const Matrix<T> &__matrix){  
  int nrows = __matrix.NumRows();
  int ncols = __matrix.NumCols();
  DiracGamma gamma;
  Matrix<typename MultType<T,complex<double> >::Type> bar(ncols,nrows);

  if(ncols == 1 && nrows == 4) bar = (__matrix.Adjoint()*gamma[0]);  
  else if(ncols == 4 && nrows == 4) 
    bar = (gamma[0]*(__matrix.Adjoint())*gamma[0]);
  else bar = (__matrix.Adjoint());

  return bar;
}
//_____________________________________________________________________________

#endif /* _DiracMatricies_H */
