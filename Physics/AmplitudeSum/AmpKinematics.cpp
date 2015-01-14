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
_motherRadius(other._motherRadius),
_wignerD(other._wignerD)
{

};

AmpKinematics::AmpKinematics(std::shared_ptr<DoubleParameter> mR, int subSys, int spin, int m, int n,
		barrierType type, std::shared_ptr<DoubleParameter> mesonRadius, std::shared_ptr<DoubleParameter> motherRadius) :
		_M(-999), _mR(mR),_subSys(subSys), _type(type), _spin(spin),_m(m),_n(n),
		_mesonRadius(mesonRadius), _motherRadius(motherRadius), _wignerD(subSys,spin)
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

void AmpKinematics::setDecayMasses(double ma, double mb, double mc, double M) {
	_ma = ma;
	_mb = mb;
	_mc = mc;
	_M = M;
}
double AmpKinematics::qValue(double x, double ma, double mb){
	double mapb = ma + mb;
	double mamb = ma - mb;
	double xSq =x*x;
	double t1 = xSq - mapb*mapb;
	double t2 = xSq - mamb*mamb;

	if( t1 < 0 ) {
		//std::cout<<"AmpKinematics: Trying to calculate break-up momentum below threshold!"<<std::endl;
		return 1; //below threshold
	}
	double result=sqrt( t1 * t2 ) / (2. * x );
	return result;

}
double AmpKinematics::BlattWeiss(double x, double mR, double ma, double mb, double spin, double mesonRadius){
	double qR = qValue(mR,ma,mb);
	double qX = qValue(x,ma,mb);
	double t0= qR*qR * mesonRadius*mesonRadius;
	double t=  qX*qX* mesonRadius*mesonRadius;
	return FormFactor(t0,t,spin);
}

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
