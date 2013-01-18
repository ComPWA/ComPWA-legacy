// PolVector class source file
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
#include "PolVector.h"
//_____________________________________________________________________________
/** @file PolVector.C
 *  @brief PolVector class source file
 */
//_____________________________________________________________________________

void PolVector::_Init(const Spin &__spin){
  _spin = __spin;
  int size = (int)(2*__spin + 1);
  _pols.resize(size);
  for(int i = 0; i < size; i++) _pols[i].SetRank(__spin);
  _projector.SetRank((int)(2*__spin));
}
//_____________________________________________________________________________

void PolVector::SetP4(const Vector4<double> &__p4,double __mass){

  this->Zero();
  _mass = __mass;
  _p4 = __p4;

  double bx = __p4.X()/__p4.E();
  double by = __p4.Y()/__p4.E();
  double bz = __p4.Z()/__p4.E();

  // indicies in pol are i = S + Mz
  if(_spin == 1) {
    _pols[2](1) = complex<double>(-1.0/sqrt(2.),0.); // m=+1  x
    _pols[2](2) = complex<double>(0.,-1.0/sqrt(2.));  // m=+1 y
    _pols[0](1) = complex<double>(1.0/sqrt(2.),0.); // m=-1 x
    _pols[0](2) = complex<double>(0.,-1.0/sqrt(2.)); // m=-1 y
    _pols[1](3) = complex<double>(1.,0.); // m=0 z

    if(__mass > 0.) this->_BoostPolVectors(bx,by,bz);    
    else{ // photon case
      _pols[1].Zero();  // mz = 0    
      if((abs(bx) > 1.E-4)||(abs(by) > 1.E-4)){
	complex<double> i(0,1);
	double x = -by/sqrt(bx*bx + by*by);
	double y = bx/sqrt(bx*bx + by*by);
	double c = _p4.CosTheta(),s = sin(_p4.Theta());
	_pols[2](0) = 0; // t
	_pols[2](1) = (-(x*x*(1-c)+c) - i*x*y*(1-c))/sqrt(2); // x
	_pols[2](2) = (-x*y*(1-c)-i*(y*y*(1-c)+c))/sqrt(2); // y
	_pols[2](3) = (y*s-i*x*s)/sqrt(2);
	_pols[0](0) = 0; // t
	_pols[0](1) = ((x*x*(1-c)+c) - i*x*y*(1-c))/sqrt(2); // x
	_pols[0](2) = (x*y*(1-c)-i*(y*y*(1-c)+c))/sqrt(2); // y
	_pols[0](3) = (-y*s-i*x*s)/sqrt(2);
      }  
    }
  }
  else{
    Spin m1,mJ,M;
    PolVector eps1(1); // spin-1 polarization vector
    eps1.SetP4(__p4,__mass);
    Spin J = _spin - 1;
    PolVector epsJ(J); // spin S-1 polarization vector
    epsJ.SetP4(__p4,__mass);

    for(m1 = -1; m1 <= 1; m1++){
      for(mJ = -J; mJ <=J; mJ++){
	for(int i = 0; i < (2*_spin + 1); i++){
	  M = i - _spin; // mz
	  _pols[i] += Clebsch(1,m1,J,mJ,_spin,M)*(eps1(m1)%epsJ(mJ));
	}
      }
    }
  }  
  this->_SetProjector(); 
}
//_____________________________________________________________________________

void PolVector::_SetProjector(){
  MetricTensor g;

  if(_spin == 1){
    if(_mass > 0.) _projector = (_p4%_p4)/(_mass*_mass) - g;    
    else _projector = -1.*g;    
  }
  else{
    _projector.Zero();
    int size = (int)(2*_spin + 1);
    for(int i = 0; i < size; i++) _projector += _pols[i]%_pols[i].Conjugate(); 
  }
}
//_____________________________________________________________________________

Tensor<complex<double> > PolVector::Propagator() const {

  Tensor<complex<double> > prop((int)(2*_spin));  
  prop = (this->_projector)/(_p4.M2() - _mass*_mass);
  return prop;
}
//_____________________________________________________________________________

Tensor<complex<double> > PolVector::PropagatorBW(double __width) const {

  Tensor<complex<double> > prop((int)(2*_spin));
  prop = (this->_projector)/(_p4.M2() - _mass*_mass 
			     + complex<double>(0.,_mass*__width));
  return prop;
}
//_____________________________________________________________________________

void PolVector::Projector(const Spin &__j,int __rank,
			  const Vector4<double> &__p4,
			  double __mass,Tensor<complex<double> > &__projector){
  
  if(__j.Denominator() == 2){
    cout << "Error! <PolVector::Projector> Spin must be integral." << endl;
    abort();
  }
  __projector.SetRank(2*__rank);
  __projector.Zero();

  Spin j = __rank - 1;
  PolVector epsj(j),eps1(1);
  eps1.SetP4(__p4,__mass);
  epsj.SetP4(__p4,__mass);

  int num_states = (int)(2*__j+1);
  for(int i = 0; i < num_states; i++) {
    Tensor<complex<double> > epsJ(__rank);
    Spin m = -__j + i;
    for(Spin m1 = -1; m1 <= 1; m1++){
      for(Spin mj = -j; mj <= j; mj++)
	epsJ += Clebsch(1,m1,j,mj,__j,m)*(eps1(m1)%epsj(mj));
    }
    __projector += epsJ%(epsJ.Conjugate());
  }
}
//_____________________________________________________________________________
