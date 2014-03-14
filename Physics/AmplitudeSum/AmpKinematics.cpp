//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

#include "Physics/AmplitudeSum/AmpKinematics.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

AmpKinematics::AmpKinematics(const AmpKinematics& other) :
_M(other._M),
_mR(other._mR),
_type(other._type),
_spin(other._spin),
_m(other._m),
_n(other._n),
_mesonRadius(other._mesonRadius),
_motherRadius(other._motherRadius)
{

};
AmpKinematics::AmpKinematics(DoubleParameter mR, int subSys, int spin, int m, int n,
		barrierType type, double mesonRadius=1.5, double motherRadius=1.5) :
		_M(-999), _mR(mR),_subSys(subSys), _type(type), _spin(spin),_m(m),_n(n),
		_mesonRadius(mesonRadius), _motherRadius(motherRadius)
{
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	_M=kin->M;
	if(_subSys==5){
		_ma=kin->m3;
		_mb=kin->m2;
		_mc=kin->m1;}
	if(_subSys==4){
		_ma=kin->m3;
		_mb=kin->m1;
		_mc=kin->m2;}
	if(_subSys==3){
		_ma=kin->m2;
		_mb=kin->m1;
		_mc=kin->m3;}
}

AmpKinematics::AmpKinematics(double ma, double mb , double mc, double M, DoubleParameter mR,
		int subSys, barrierType type, int spin, int m, int n, double mesonR=1.5, double motherR=1.5) :
				_ma(ma), _mb(mb), _mc(mc), _M(M), _mR(mR),_subSys(subSys), _type(type),
				_spin(spin),_m(m),_n(n), _mesonRadius(mesonR), _motherRadius(motherR)
{

};

void AmpKinematics::setDecayMasses(double ma, double mb, double mc, double M) {
	_ma = ma;
	_mb = mb;
	_mc = mc;
	_M = M;
}
double AmpKinematics::qValue(double x, double ma, double mb){
	double mapb = ma + mb;
	double mamb = ma - mb;

	if( (x*x - mapb*mapb) < 0 ) {
		//std::cout<<"AmpKinematics: Trying to calculate break-up momentum below threshold!"<<std::endl;
		return 1; //below threshold
	}
	double result=sqrt( (x*x - mapb*mapb) * (x*x - mamb*mamb) ) / (2. * x );
	//	std::cout<<"ss "<<ma<<" "<<mb<<" "<<x<<" "<<sqrt(x*x/4-ma*ma)<<" " <<result<<std::endl;
	return result;

}
//double AmpKinematics::q0(double ma, double mb) const {
//	double mapb = ma + mb;
//	double mamb = ma - mb;
//
//	double mr = _mR.GetValue();
//
//	if( (mr*mr - mapb*mapb) < 0 ) {
//		//std::cout<<"AmpKinematics: Trying to calculate break-up momentum below threshold!"<<std::endl;
//		return 1; //below threshold
//	}
//	return sqrt( (mr*mr - mapb*mapb) * (mr*mr - mamb*mamb)) / (2. * mr );
//}
//double AmpKinematics::q(double x, double ma, double mb) const {
//	double mapb = ma + mb;
//	double mamb = ma - mb;
//
//	if( (x*x - mapb*mapb) < 0 ) {
//		//std::cout<<"AmpKinematics: Trying to calculate break-up momentum below threshold!"<<std::endl;
//		return 1; //below threshold
//	}
//	double result=sqrt( (x*x - mapb*mapb) * (x*x - mamb*mamb) ) / (2. * x );
//	//	std::cout<<"ss "<<ma<<" "<<mb<<" "<<x<<" "<<sqrt(x*x/4-ma*ma)<<" " <<result<<std::endl;
//	return result;
//}

// compute square of Blatt-Weisskopf barrier factor
//double AmpKinematics::BLres2(double x) const {
//
//	double t0= q0()*q0() * _mesonRadius*_mesonRadius;
//	double t= q(x)*q(x) * _mesonRadius*_mesonRadius;
//	return FormFactor(t0,t);
//}
double AmpKinematics::BlattWeiss(double x, double mR, double ma, double mb, double spin, double mesonRadius){
	double t0= qValue(mR,ma,mb)*qValue(mR,ma,mb) * mesonRadius*mesonRadius;
	double t= qValue(x,ma,mb)*qValue(x,ma,mb) * mesonRadius*mesonRadius;
	return FormFactor(t0,t,spin);
}
//double AmpKinematics::BLmother2(double x) const {
//
//	//calculate q* here
//	double t0= q0()*q0() * _motherRadius*_motherRadius;
//	double t= q(x)*q(x) * _motherRadius*_motherRadius;
//	return FormFactor(t0,t);
//}
double AmpKinematics::FormFactor(double z0, double z,unsigned int spin){
	double nom=0, denom=0;
	if (spin == 0) return 1;
	else if (spin == 1){
		//		if(_type==barrierType::BWPrime){
		nom = 1 + z0;
		denom = 1 + z;
		//		} else if(_type==barrierType::BW){
		//			nom = 2*z;
		//			denom = 1 + z;
		//		} else {
		//			std::cout<<"Wrong BLW factor definition: "<<_type<<std::endl;
		//			return 1;
		//		}
	}
	else if (spin == 2) {
		//		if(_type==barrierType::BWPrime){
		nom = (z0-3)*(z0-3)+9*z0;
		denom = (z-3)*(z-3)+9*z;
		//		} else if(_type==barrierType::BW){
		//			nom = 13*z*z;
		//			denom = (z-3)*(z-3)+9*z;
		//		} else {
		//			std::cout<<"Wrong BLW factor definition: "<<_type<<std::endl;
		//			return 1;
		//		}
	}
	else{
		std::cout<<"Wrong spin value! BLW factors only implemented for spin 0,1 and 2! "<<std::endl;
	}
	return nom/denom;
}
