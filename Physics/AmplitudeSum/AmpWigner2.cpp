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
	DalitzKinematics* kin = DalitzKinematics::instance();
	_M=kin->M;
	_m1=kin->m1;
	_m2=kin->m2;
	_m3=kin->m3;

	_spinM = kin->getSpin(0);
	_spin1 = kin->getSpin(1);
	_spin2 = kin->getSpin(2);
	_spin3 = kin->getSpin(3);

}
double AmpWigner2::evaluate() const {
	dataPoint* point = dataPoint::instance();
	DalitzKinematics* kin = DalitzKinematics::instance();

	double cosTheta=-999, result=-999;
	Spin J((int)_resSpin);
	/*
	 * \in and \out are the difference in the spins of the ingoing and outgoing particles.
	 * In literature it is often denoted with mu and muPrime
	 */
	Spin out, in;

	switch(_subSys){
	case 3:
		cosTheta = kin->calcHelicityAngle(point->getMsq(3),point->getMsq(4),_M,_m3,_m1,_m2);
		in = (int)(_spinM-_spin1); out = (int)(_spin3-_spin2);
		break;
	case 4:
		cosTheta = kin->calcHelicityAngle(point->getMsq(4),point->getMsq(5),_M,_m2,_m1,_m3);
		in = (int)(_spinM-_spin2); out = (int)(_spin3-_spin1);
		break;
	case 5:
		cosTheta = kin->calcHelicityAngle(point->getMsq(5),point->getMsq(4),_M,_m1,_m2,_m3);
		in = (int)(_spinM-_spin3); out = (int)(_spin1-_spin2);
		break;
	default:
		std::cout<<"AmpWigner: wrong subSystem! Exit!"<<std::endl; exit(1);
	}
	double theta = acos(cosTheta);

	/*
	 * Calling WignerD function
	 * Note that Wigner_d depends on the sign of \in and \out. I hope it is correctly assigned.
	 */
	result = Wigner_d(J,in,out,theta);
	if( ( result!=result ) || (theta!=theta)) {
		std::cout<< "NAN! J="<< J<<" M="<<out<<" N="<<in<<" beta="<<theta<<std::endl;
		return 0;
	}
	return result;
}
