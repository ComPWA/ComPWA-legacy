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
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
AmpFlatteRes::AmpFlatteRes(const char *name,
		std::shared_ptr<DoubleParameter> resMass,
		std::shared_ptr<DoubleParameter> mesonRadius, //  meson radius
		std::shared_ptr<DoubleParameter> motherRadius, //  mother radius
		std::shared_ptr<DoubleParameter> g1, std::shared_ptr<DoubleParameter> g2,
		double g2_partA, double g2_partB,
		int subSys, ///  meson radius
		int resSpin, int m, int n) :
		AmpAbsDynamicalFunction(name),
		AmpKinematics(resMass, subSys, resSpin, m, n, AmpKinematics::barrierType(BWPrime),
				mesonRadius, motherRadius),
		_g2(g2), _g1(g1),
		_g2_partA(g2_partA), _g2_partB(g2_partB),
		foundMasses(false),
		_wignerD(subSys,resSpin),
		nParams(5)
{
}

AmpFlatteRes::AmpFlatteRes(const AmpFlatteRes& other, const char* newname) :
						  AmpAbsDynamicalFunction(other, newname),
						  AmpKinematics(other),
						  _g2(other._g2),
						  _g1(other._g1),
						  _wignerD(other._wignerD)
{
}

AmpFlatteRes::AmpFlatteRes(const AmpFlatteRes& other) :
						  AmpAbsDynamicalFunction(other),
						  AmpKinematics(other),
						  _g2(other._g2),
						  _g1(other._g1),
						  _wignerD(other._wignerD)
{
}

AmpFlatteRes::~AmpFlatteRes() 
{
}

void AmpFlatteRes::setDecayMasses(double ma, double mb, double mc, double M){
	AmpKinematics::setDecayMasses(ma, mb, mc, M);
	return;
}

void AmpFlatteRes::setBarrierMass(double mBarA, double mBarB) {
	_g2_partA = mBarA;
	_g2_partB = mBarB;
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
std::complex<double> AmpFlatteRes::dynamicalFunction(double mSq, double mR, double ma, double mb, double g1,
		double mHiddenA, double mHiddenB, double g2,unsigned int J ){
	double p1 = 2*AmpKinematics::qValue(sqrt(mSq), ma,mb)/sqrt(mSq);//break-up momenta decay channel (e.g. a0->KK)
	double p2 = 2*AmpKinematics::qValue(sqrt(mSq), mHiddenA,mHiddenB)/sqrt(mSq);//break-up momenta hidden channel (e.g. a0->eta pi)

	std::complex<double> denom(mR*mR - mSq, -p1*g1*g1-p2*g2*g2);

	std::complex<double> result = (std::complex<double>( g1*g1 , 0 ) / denom);

	if(result.real()!=result.real()) {std::cout << "RE part NAN" << std::endl; return 0;}
	if(result.imag()!=result.imag()) {std::cout << "IM part NAN" << std::endl; return 0;}
	return result;
}
