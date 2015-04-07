//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding correct g1s
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include <math.h>
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
AmpFlatteRes::AmpFlatteRes(const char *name,
		std::shared_ptr<DoubleParameter> resMass,
		std::shared_ptr<DoubleParameter> mesonRadius, //  meson radius
		std::shared_ptr<DoubleParameter> motherRadius, //  mother radius
		std::shared_ptr<DoubleParameter> g1, std::shared_ptr<DoubleParameter> g2,
		double g2_partA, double g2_partB, int subSys, int resSpin, int m, int n) :
		AmpAbsDynamicalFunction(name),
		AmpKinematics(resMass, subSys, resSpin, m, n, mesonRadius, motherRadius),
		_g2(g2), _g1(g1), _g2_partA(g2_partA), _g2_partB(g2_partB),
		foundMasses(false),	nParams(5)
{
}

AmpFlatteRes::~AmpFlatteRes() 
{
}

std::complex<double> AmpFlatteRes::evaluateAmp(dataPoint& point) {
	if(_g2_partA<0||_g2_partA>5||_g2_partB<0||_g2_partB>5) {
		cout<<"Barrier masses not set! Use setBarrierMass() first!"<<endl;
		return 0;
	}
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	if(!foundMasses){
		id23 = point.getID("m23sq");
		id13 = point.getID("m13sq");
		foundMasses=true;
	}
	double mSq=-999;
	switch(_subSys){
	case 3:
		mSq=(kin->getThirdVariableSq(point.getVal(0),point.getVal(1)));	break;
	case 4: mSq=(point.getVal(1)); break;
	case 5: mSq=(point.getVal(0)); break;
	}

	return dynamicalFunction(mSq,_mR->GetValue(),_ma,_mb,_g1->GetValue(), _g2_partA,_g2_partB,_g2->GetValue(),_spin);
}
std::complex<double> AmpFlatteRes::dynamicalFunction(double mSq, double mR,
		double massA1, double massA2, double gA,
		double massB1, double massB2, double gB,
		unsigned int J ){
	double sqrtS = sqrt(mSq);
	double den = 8*M_PI*sqrtS; //use Pi value of math.h

	double mesonRadius = 1.5;//TODO: pass this as argument

	//channel A - signal channel
	//break-up momentum
	double pA = AmpKinematics::qValue(sqrtS, massA1,massA2) / den;
	double barrierA = AmpKinematics::FormFactor(sqrtS,mR,massA1,massA2,J,mesonRadius)/AmpKinematics::FormFactor(mR,mR,massA1,massA2,J,mesonRadius);
	double qTermA = std::pow((qValue(sqrtS,massA1,massA2) / qValue(mR,massA1,massA2)), (2.*J+ 1.));
	//convert coupling to partial width of channel A
	double gammaA = couplingToWidth(mSq,mR,gA,massA1,massA2,J,mesonRadius);
	double termA = gammaA*qTermA*barrierA*barrierA;
	//channel B - hidden channel
	//break-up momentum
	double pB = AmpKinematics::qValue(sqrtS, massB1,massB2) / den;
	double barrierB = AmpKinematics::FormFactor(sqrtS,mR,massB1,massB2,J,1.5)/AmpKinematics::FormFactor(mR,mR,massB1,massB2,J,1.5);
	double qTermB = std::pow((qValue(sqrtS,massB1,massB2) / qValue(mR,massB1,massB2)), (2.*J+ 1.));
	//convert coupling to partial width of channel B
	double gammaB = couplingToWidth(mSq,mR,gB,massB1,massB2,J,mesonRadius);
	double termB = gammaB*qTermB*barrierB*barrierB;

	//Coupling constant from production reaction. In case of a particle decay the production
	//coupling doesn't depend in energy since the CM energy is in the (RC) system fixed to the
	//mass of the decaying particle
	double g_production = 1;

	//old approach
//	std::complex<double> denom( mR*mR - mSq, -pA*gA*gA-pB*gB*gB );
	//new approach
	std::complex<double> denom( mR*mR - mSq, (-1)*sqrtS*(termA + termB) );
	//std::cout<<denom<<" "<<denom2<<" "<<barrierA<< " "<<barrierB<<std::endl;
	std::complex<double> result = (std::complex<double>( gA*g_production , 0 ) / denom);

	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"nan in Flatte: "<<barrierA<<" "<<mR<<" "<<mSq<<" "<<massA1<<" "<<massA2<<std::endl;
		return 0;
	}
	return result;
}
