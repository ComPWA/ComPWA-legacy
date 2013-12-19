//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff  -
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

#include <cmath>
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

#include "qft++.h"

AmpWigner2::AmpWigner2(unsigned int subSys, unsigned int resSpin) : _resSpin(resSpin),_subSys(subSys)
{
	initialise();
}

AmpWigner2::AmpWigner2(const AmpWigner2& other, const char* newname) : _resSpin(other._resSpin),_subSys(other._subSys)
{
	initialise();
}

void AmpWigner2::initialise()
{
	static dataPoint* point = dataPoint::instance();
	_M=point->DPKin.M;
	_m1=point->DPKin.m1;
	_m2=point->DPKin.m2;
	_m3=point->DPKin.m3;

	_spinM = point->DPKin.getSpin(0);
	_spin1 = point->DPKin.getSpin(1);
	_spin2 = point->DPKin.getSpin(2);
	_spin3 = point->DPKin.getSpin(3);

}
double AmpWigner2::evaluate() const {
	dataPoint* point = dataPoint::instance();
	DPKinematics kin = point->DPKin;

	double cosTheta=-999, result=-999;
	Spin J((int)_resSpin);
	/*
	 * \in and \out are the difference in the spins of the ingoing and outgoing particles.
	 * In literature it is often denoted with mu and muPrime
	 */
	Spin out, in;

	switch(_subSys){
	case 3:
		cosTheta = kin.calcHelicityAngle(point->getMsq(3),point->getMsq(4),_M,_m3,_m1,_m2);
		in = (int)(_spinM-_spin3); out = (int)(_spin1-_spin2);
		break;
	case 4:
		cosTheta = kin.calcHelicityAngle(point->getMsq(4),point->getMsq(5),_M,_m2,_m3,_m1);
		in = (int)(_spinM-_spin2); out = (int)(_spin3-_spin1);
		break;
	case 5:
		cosTheta = kin.calcHelicityAngle(point->getMsq(5),point->getMsq(4),_M,_m1,_m2,_m3);
		in = (int)(_spinM-_spin1); out = (int)(_spin3-_spin2);
		break;
	default:
		std::cout<<"AmpWigner: wrong subSystem! Exit!"<<std::endl; exit(1);
	}
	if(cosTheta>1.) cosTheta=1.;
	if(cosTheta<-1.) cosTheta=-1.;
	double theta = acos(cosTheta);

	/*
	 * Calling WignerD function
	 * Note that Wigner_d depends on the sign of \in and \out. I hope it is correctly assigned.
	 */
	result = Wigner_d(J,in,out,theta);
	if( ( result!=result ) || (theta!=theta)) {
		std::cout<< "NAN! J="<< J<<" M="<<out<<" N="<<in<<" beta="<<cosTheta<<std::endl;
		std::cout<< "msq12="<< point->getMsq(3)<<" msq13="<< point->getMsq(4)<<" msq23="<< point->getMsq(5)<<std::endl;
		return 0;
	}
	return result;
}

/*double AmpWigner2::evaluateTree(const ParameterList& paras, const std::string name) const{
	dataPoint* point = dataPoint::instance();
	DPKinematics kin = point->DPKin;

	//double Gamma0, GammaV;
	double m23 = double(paras.GetParameterValue("m23"));
	double m13 = double(paras.GetParameterValue("m13"));
	double m12 = double(paras.GetParameterValue("m12"));
	double M  = double(paras.GetParameterValue("M"));
	double m1 = double(paras.GetParameterValue("m1"));
	double m2 = double(paras.GetParameterValue("m2"));
	double m3 = double(paras.GetParameterValue("m3"));
	//double locmax_sq = Double_t(paras.GetParameterValue("mb_"+name));
	//double locmin_sq = Double_t(paras.GetParameterValue("mb_"+name));
	unsigned int subSysFlag = double(paras.GetParameterValue("subSysFlag_"+name));
	double motherSpin = double(paras.GetParameterValue("motherSpin_"+name));
	double resSpin = double(paras.GetParameterValue("resSpin_"+name));
	double outSpin1 = double(paras.GetParameterValue("outSpin1_"+name));
	double outSpin2 = double(paras.GetParameterValue("outSpin2_"+name));
	double outSpin3 = double(paras.GetParameterValue("outSpin3_"+name));

	double cosTheta=-999, result=-999;
	Spin J((int)resSpin);
	Spin out, in;

	switch(_subSys){
	case 3:
		cosTheta = kin.calcHelicityAngle(m12,m13,M,m3,m1,m2);
		in = (int)(motherSpin-outSpin3); out = (int)(outSpin1-outSpin2);
		break;
	case 4:
		cosTheta = kin.calcHelicityAngle(m13,m23,M,m2,m3,m1);
		in = (int)(motherSpin-outSpin2); out = (int)(outSpin3-outSpin1);
		break;
	case 5:
		cosTheta = kin.calcHelicityAngle(m23,m13,M,m1,m2,m3);
		in = (int)(motherSpin-outSpin1); out = (int)(outSpin3-outSpin2);
		break;
	default:
		std::cout<<"AmpWigner: wrong subSystem! Exit!"<<std::endl; exit(1);
	}
	if(cosTheta>1.) cosTheta=1.;
	if(cosTheta<-1.) cosTheta=-1.;
	double theta = acos(cosTheta);

	result = Wigner_d(J,in,out,theta); //TODO: use same functions as above
	if( ( result!=result ) || (theta!=theta)) {
		std::cout<< "NAN! J="<< J<<" M="<<out<<" N="<<in<<" beta="<<cosTheta<<std::endl;
		std::cout<< "msq12="<< m12<<" msq13="<< m13<<" msq23="<< m23<<std::endl;
		return 0;
	}
	return result;
}*/
