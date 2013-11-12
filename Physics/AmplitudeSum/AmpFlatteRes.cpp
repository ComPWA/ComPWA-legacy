//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding correct couplings
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
//#include "RooRealVar.h"
AmpFlatteRes::AmpFlatteRes(const char *name,
		DoubleParameter& resMass, DoubleParameter& resWidth,
		double& mesonRadius, ///  meson radius
		DoubleParameter& couplingHidden, DoubleParameter& coupling,
		double massHiddenChannelA, double massHiddenChannelB,
		int subSys, ///  meson radius
		int resSpin, int m, int n) :
		AmpAbsDynamicalFunction(name),
		AmpKinematics(resMass, subSys, resSpin, m, n, AmpKinematics::barrierType(BWPrime), mesonRadius, 1.5),
		_couplingHiddenChannel(couplingHidden),
		_coupling(coupling),
		_massHiddenChannelA(massHiddenChannelA),
		_massHiddenChannelB(massHiddenChannelB),
		_wignerD(name, resSpin,m,n, subSys)
{
	initialise();
}

AmpFlatteRes::AmpFlatteRes(const AmpFlatteRes& other, const char* newname) :
		  AmpAbsDynamicalFunction(other, newname),
		  AmpKinematics(other),
		  _couplingHiddenChannel(other._couplingHiddenChannel),
		  _coupling(other._coupling),
		  _wignerD(other._wignerD)
{
	initialise();
}

AmpFlatteRes::AmpFlatteRes(const AmpFlatteRes& other) :
		  AmpAbsDynamicalFunction(other),
		  AmpKinematics(other),
		  _couplingHiddenChannel(other._couplingHiddenChannel),
		  _coupling(other._coupling),
		  _wignerD(other._wignerD)
{
	initialise();
}

AmpFlatteRes::~AmpFlatteRes() 
{
}

void AmpFlatteRes::initialise() 
{
}
  void AmpFlatteRes::setDecayMasses(double ma, double mb, double mc, double M){
	  AmpKinematics::setDecayMasses(ma, mb, mc, M);
	  _wignerD.setDecayMasses(ma, mb, mc, M);
	  return;
  }

void AmpFlatteRes::setBarrierMass(double mBarA, double mBarB) {
	_massHiddenChannelA = mBarA;
	_massHiddenChannelB = mBarB;
}

std::complex<double> AmpFlatteRes::evaluate() const {
	if(_massHiddenChannelA<0||_massHiddenChannelA>5||_massHiddenChannelB<0||_massHiddenChannelB>5) {
		cout<<"Barrier masses not set! Use setBarrierMass() first!"<<endl;
		return 0;
	}
	//	double m0 = Double_t(_m0);
	double m = dataPoint::instance()->getM(_subSys);
	double spinTerm = _wignerD.evaluate();
//	double m  = Double_t(_x23);

	double p1 = 2*q(m, _massHiddenChannelA,_massHiddenChannelB)/m;//break-up momenta hidden channel (e.g. a0->eta pi)
	double p2 = 2*q(m, _ma,_mb)/m;//break-up momenta decay channel (e.g. a0->KK)
	double g1 = _couplingHiddenChannel.GetValue();//couppling a0->eta pi
	double g2 = _coupling.GetValue();//coupling a0->KK

	std::complex<double> denom(_mR*_mR - m*m, -p1*g1*g1-p2*g2*g2);

	//	RooComplex result = (RooComplex(g2*g2,0) / denom); //use KK decay channel here
	std::complex<double> result = (std::complex<double>(_norm * spinTerm * g2*g2,0) / denom); //use KK decay channel here

	if(result.real()!=result.real()) {std::cout << "RE part NAN" << std::endl; return 0;}
	if(result.imag()!=result.imag()) {std::cout << "IM part NAN" << std::endl; return 0;}
	return result;

}
