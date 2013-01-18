// DiracAntiSpinor class source file
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
#include "DiracAntiSpinor.h"
//_____________________________________________________________________________
/** @file DiracAntiSpinor.C
 *  @brief DiracAntiSpinor class source file
 */
//_____________________________________________________________________________

void DiracAntiSpinor::SetP4(const Vector4<double> &__p4,double __mass){

  Matrix <complex<double> > sigP(2,2);
  Matrix <complex<double> > chi(2,1);
  PauliSigma sigma;

  _p4 = __p4;
  _mass = __mass;

  complex<double> norm = sqrt(__p4.E() + __mass);
  complex<double> epm = __p4.E() + __mass;
  sigP = sigma[1]*__p4.X() + sigma[2]*__p4.Y() + sigma[3]*__p4.Z();

  // spin up
  chi(0,0) = 0.;  
  chi(1,0) = 1.;  
  _spinors[1](2,0) = chi(0,0)*norm;
  _spinors[1](3,0) = chi(1,0)*norm;
  _spinors[1](0,0) = ((sigP*chi)(0,0))*norm/epm;
  _spinors[1](1,0) = ((sigP*chi)(1,0))*norm/epm;

  // spin down
  chi(0,0) = 1.;
  chi(1,0) = 0.;
  _spinors[0](2,0) = chi(0,0)*norm;
  _spinors[0](3,0) = chi(1,0)*norm;
  _spinors[0](0,0) = ((sigP*chi)(0,0))*norm/epm;
  _spinors[0](1,0) = ((sigP*chi)(1,0))*norm/epm;
  
  this->_SetProjector();
}
//_____________________________________________________________________________

void DiracAntiSpinor::Boost(double __bx,double __by,double __bz) {

  _p4.Boost(__bx,__by,__bz);
  this->_SetProjector();

  double beta = sqrt(__bx*__bx + __by*__by + __bz*__bz);
  double ux = __bx/beta;
  double uy = __by/beta;
  double uz = __bz/beta;
  double gamma = 1.0/sqrt(1. - beta*beta);
  double th = sqrt((gamma - 1.0)/(gamma + 1.0));
  double ch = 1.0/sqrt(1. - th*th);
  double sh = ch*th;

  // get the spin 1/2 boost operator
  Matrix<complex<double> > D(4,4);

  D(0,0) = ch;
  D(1,1) = ch;
  D(2,2) = ch;
  D(3,3) = ch;
  D(0,2) = uz*sh;
  D(0,3) = (ux + complex<double>(0.,-uy))*sh;
  D(1,2) = (ux + complex<double>(0.,uy))*sh;
  D(1,3) = -uz*sh;
  D(2,0) = D(0,2);
  D(2,1) = D(0,3);
  D(3,0) = D(1,2);
  D(3,1) = D(1,3);
  
  for(int i = 0; i < 2; i++) _spinors[i] = D*_spinors[i];
}
//_____________________________________________________________________________
