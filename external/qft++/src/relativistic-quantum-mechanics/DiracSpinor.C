// DiracSpinor class source file
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
#include "DiracSpinor.h"
//_____________________________________________________________________________
/** @file DiracSpinor.C
 *  @brief DiracSpinor class source file
 */
//_____________________________________________________________________________

void DiracSpinor::_Init(const Spin &__spin) {
  _spin = __spin;
  int two_S = (int)(2*__spin);
  _spinors.resize(two_S + 1);
  _projector.Resize(4,4);
  SetRank(_projector,(two_S - 1));
  for(int i = 0; i < (int)_spinors.size(); i++){
    _spinors[i].Resize(4,1);
    SetRank(_spinors[i],(int)((two_S - 1)/2));
  }
}
//_____________________________________________________________________________

void DiracSpinor::_Copy(const DiracSpinor &__dspin){
  _spin = __dspin._spin;
  _mass = __dspin._mass;
  _p4 = __dspin._p4;
  int size = (int)__dspin._spinors.size();
  this->_Init(__dspin._spin);
  _projector = __dspin._projector;
  for(int i = 0; i < size; i++) _spinors[i] = __dspin._spinors[i];
}
//_____________________________________________________________________________

void DiracSpinor::SetP4(const Vector4<double> &__p4,double __mass){

  this->Zero();
  _p4 = __p4;
  _mass = __mass;
  if(_spin == 1/2.){ // spin-1/2 case 
    complex<double> norm,epm;
    Matrix<complex<double> > sigP(2,2);
    Matrix<complex<double> > chi(2,1);
    PauliSigma sigma;

    norm = sqrt(__p4.E() + __mass);
    epm = __p4.E() + __mass;    
    sigP = sigma(1)*__p4.X() + sigma(2)*__p4.Y() + sigma(3)*__p4.Z();

    // set spin up    
    chi(0,0) = 1.;  
    chi(1,0) = 0.;
    _spinors[1](0,0) = chi(0,0)*norm;
    _spinors[1](1,0) = chi(1,0)*norm;
    _spinors[1](2,0) = ((sigP*chi)(0,0))*(norm/epm);
    _spinors[1](3,0) = ((sigP*chi)(1,0))*(norm/epm);

    // set spin down    
    chi(0,0) = 0.;
    chi(1,0) = 1.;
    _spinors[0](0,0) = chi(0,0)*norm;
    _spinors[0](1,0) = chi(1,0)*norm;
    _spinors[0](2,0) = ((sigP*chi)(0,0))*(norm/epm);
    _spinors[0](3,0) = ((sigP*chi)(1,0))*(norm/epm);
  }
  else { // higher spins
    Spin j = _spin - 1/2.;
    DiracSpinor u; // spin-1/2
    PolVector eps(j); 
    u.SetP4(__p4,__mass);
    eps.SetP4(__p4,__mass);
    for(Spin mf = -1/2.; mf <= 1/2.; mf++){
      for(Spin mb = -j; mb <= j; mb++){
	for(int i = 0; i < (int)(2*_spin + 1); i++){
	  Spin mz = i - _spin;
	  _spinors[i] += u(mf)*Clebsch(j,mb,1/2.,mf,_spin,mz)*eps(mb);
	}
      }
    }
  }
  this->_SetProjector();
}
//_____________________________________________________________________________

void DiracSpinor::_SetProjector() {
  if(_spin == 1/2.){
    IdentityMatrix<double> I(4);
    DiracGamma gamma;    
    _projector = (_p4*gamma + I*_mass)/(2.*_mass); // ("p slash" + m)/2m
  }
  else{ // higher spin    
    _projector.Zero();
    int size = (int)(2*_spin + 1);
    for(int i = 0; i < size; i++) _projector += _spinors[i] % Bar(_spinors[i]);
    _projector /= (2*_mass);
  }
}
//_____________________________________________________________________________
  
void DiracSpinor::Boost(double __bx,double __by,double __bz){

  if(_spin == 1/2.){ 
    double beta = sqrt(__bx*__bx + __by*__by + __bz*__bz);
    double ux = __bx/beta;
    double uy = __by/beta;
    double uz = __bz/beta;
    double gamma = 1.0/sqrt(1. - beta*beta);
    double th = sqrt((gamma - 1.0)/(gamma + 1.0));
    double ch = 1.0/sqrt(1. - th*th);
    double sh = ch*th;

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
  
    _spinors[1] = D*_spinors[1];
    _spinors[0] = D*_spinors[0];
  }
  else { // higher spin ...
    this->Zero();
    Spin j = _spin - 1/2.;
    DiracSpinor u; // spin-1/2
    PolVector eps(j); 
    u.SetP4(_p4,_mass); u.Boost(__bx,__by,__bz);
    eps.SetP4(_p4,_mass); eps.Boost(__bx,__by,__bz);
    for(Spin mf = -1/2.; mf <= 1/2.; mf++){
      for(Spin mb = -j; mb <= j; mb++){
	for(int i = 0; i < (int)(2*_spin + 1); i++){
	  Spin mz = i - _spin;
	  _spinors[i] += u(mf)*Clebsch(j,mb,1/2.,mf,_spin,mz)*eps(mb);
	}
      }
    }
  }
  _p4.Boost(__bx,__by,__bz);
  this->_SetProjector();
}
//_____________________________________________________________________________

void DiracSpinor::Projector(const Spin &__j,int __rank,
			    const Vector4<double> &__p4,double __mass,
			    Matrix<Tensor<complex<double> > > &__projector){
  if(__j.Denominator() == 1){
    cout << "Error! <DiracSpinor::Projector> Spin must be half-integral." 
	 << endl;
    abort();
  }
  __projector.Resize(4,4);
  SetRank(__projector,2*__rank);
  __projector.Zero();

  Spin j = __rank;
  PolVector eps(j);
  DiracSpinor u;

  u.SetP4(__p4,__mass);
  eps.SetP4(__p4,__mass);

  int num_states = (int)(2*__j+1);
  Matrix<Tensor<complex<double> > > psi(4,1);
  SetRank(psi,__rank);

  for(int i = 0; i < num_states; i++) {
    psi.Zero();
    Spin m = -__j + i;
    for(Spin m1 = -1/2.; m1 <= 1/2.; m1++){
      for(Spin mj = -j; mj <= j; mj++){
	psi += u(m1)*Clebsch(1/2.,m1,j,mj,__j,m)*eps(mj);
      }
    }
    __projector += psi % Bar(psi);
  }
  __projector /= 2*__mass;
}
//_____________________________________________________________________________
