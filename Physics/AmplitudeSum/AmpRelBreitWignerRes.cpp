//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include <stdlib.h>


AmpRelBreitWignerRes::AmpRelBreitWignerRes(const char *name,
		std::shared_ptr<DoubleParameter> resMass, std::shared_ptr<DoubleParameter> resWidth,
		std::shared_ptr<DoubleParameter> mesonRadius, std::shared_ptr<DoubleParameter> motherRadius,
		int subSys,
		int resSpin, int m, int n) :
		AmpAbsDynamicalFunction(name),
		AmpKinematics(resMass, subSys, resSpin, m, n, mesonRadius, motherRadius),
		_resWidth(resWidth),
		foundMasses(false),
		nParams(6)
{
	initialise();
}

AmpRelBreitWignerRes::~AmpRelBreitWignerRes() 
{
}

void AmpRelBreitWignerRes::initialise() 
{
}

std::complex<double> AmpRelBreitWignerRes::evaluateAmp(dataPoint& point) {
	if(_ma == -999 || _mb == -999 ||_mc == -999 ||_M == -999){
		std::cout<<"Masses of decay products not set!"<<std::endl;
		return 0;
	}
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	if(!foundMasses){
		id23 = point.getID("m23sq");
		id13 = point.getID("m13sq");
		foundMasses=true;
	}
	double mSq = -999;
	switch(_subSys){
	case 3: mSq=kin->getThirdVariableSq(point.getVal(id23),point.getVal(id13)); break;
	case 4: mSq=point.getVal(id13); break;
	case 5: mSq=point.getVal(id23); break;
	}

	return dynamicalFunction(mSq,_mR->GetValue(),_ma,_mb,_resWidth->GetValue(),_spin,_mesonRadius->GetValue());
}
std::complex<double> AmpRelBreitWignerRes::dynamicalFunction(double mSq, double mR, double ma, double mb, double width, unsigned int J, double mesonRadius){
	double sqrtS = sqrt(mSq);

	double barrier = AmpKinematics::FormFactor(sqrtS,mR,ma,mb,J,mesonRadius)/AmpKinematics::FormFactor(mR,mR,ma,mb,J,mesonRadius);
	std::complex<double> qTerm = std::pow((qValue(sqrtS,ma,mb) / qValue(mR,ma,mb)), (2.*J+ 1.));
	//Calculate coupling constant to final state
	double g_final = widthToCoupling(mSq,mR,width,ma,mb,J,mesonRadius);

	//Coupling constant from production reaction. In case of a particle decay the production
	//coupling doesn't depend in energy since the CM energy is in the (RC) system fixed to the
	//mass of the decaying particle
	double g_production = 1;

	//	std::complex<double> denom(mR*mR - mSq, (-1)*sqrtS*(width*qTerm*barrier*barrier));
	std::complex<double> denom(mR*mR - mSq + sqrtS*(width*qTerm.imag()*barrier*barrier), (-1)*sqrtS*(width*qTerm.real()*barrier*barrier) );
	std::complex<double> result = std::complex<double>(g_final*g_production,0) / denom;

	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"nan in BW: "<<barrier<<" "<<mR<<" "<<mSq<<" "<<ma<<" "<<mb<<std::endl;
		return 0;
	}
	return result;
}
