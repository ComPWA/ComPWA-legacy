//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
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
//#include <stdlib.h>
//#include "gsl/gsl_math.h"
//#include "gsl/gsl_monte.h"
//#include "gsl/gsl_monte_plain.h"
//#include "gsl/gsl_monte_miser.h"
//#include "gsl/gsl_monte_vegas.h"

AmpRelBreitWignerRes::AmpRelBreitWignerRes(const char *name,
		DoubleParameter& resMass, DoubleParameter& resWidth,
		double& mesonRadius, //  meson radius
		int subSys,
		int resSpin, int m, int n
) :
AmpAbsDynamicalFunction(name),
AmpKinematics(resMass, subSys, resSpin, m, n, AmpKinematics::barrierType(BWPrime), mesonRadius, 1.5),
_resWidth(resWidth),
_wignerD(name, resSpin, m, n, subSys)
{
	initialise();
}

AmpRelBreitWignerRes::AmpRelBreitWignerRes(const AmpRelBreitWignerRes& other, const char* newname) :
				  AmpAbsDynamicalFunction(other, newname),
				  AmpKinematics(other),
				  _resWidth(other._resWidth),
				  _wignerD(other._wignerD)
{
	initialise();
}

AmpRelBreitWignerRes::AmpRelBreitWignerRes(const AmpRelBreitWignerRes& other) :
				  AmpAbsDynamicalFunction(other),
				  AmpKinematics(other),
				  _resWidth(other._resWidth),
				  _wignerD(other._wignerD)
{
	initialise();
}

AmpRelBreitWignerRes::~AmpRelBreitWignerRes() 
{
}

void AmpRelBreitWignerRes::initialise() 
{
}

void AmpRelBreitWignerRes::setDecayMasses(double ma, double mb, double mc, double M){
	AmpKinematics::setDecayMasses(ma, mb, mc, M);
	_wignerD.setDecayMasses(ma, mb, mc, M);
	return;
}
std::complex<double> AmpRelBreitWignerRes::evaluate() const {
	if(_ma == -999 || _mb == -999 ||_mc == -999 ||_M == -999){
		std::cout<<"Masses of decay products not set!"<<std::endl;
		return 0;
	}
	std::complex<double> result;
	double m = dataPoint::instance()->getM(_subSys);
	double spinTerm = _wignerD.evaluate();
	double Gamma0 = _resWidth.GetValue();
	double GammaV = Gamma0 * pow(q(m) / q0(), 2.*_spin + 1.) * (_mR / m) * BLres2(m);

	std::complex<double> denom(_mR*_mR - m*m, -_mR * GammaV);

	//	  result = RooComplex(spinTerm*_mR * Gamma0) / denom; //wrong!
	//  result = RooComplex(spinTerm*BLprime2(m)) / denom;
	result = std::complex<double>( spinTerm ) / denom; //Laura++ (old) definition - is this used in DKsKK analysis?
	//	result = RooComplex( spinTerm*sqrt(BLres2(m))*sqrt(BLmother2(m)) ) / denom; //Laura++ (new) definition

	if(result.real()!=result.real()) {std::cout << "RE part NAN" << std::endl;return 0;}
	if(result.imag()!=result.imag()) {std::cout << "IM part NAN" << std::endl; return 0;}
	return result;
}
double AmpRelBreitWignerRes::eval(double x[], int dim, void * param) const {

	//set values here
	dataPoint::instance()->setMsq(4,x[0]);
	dataPoint::instance()->setMsq(5,x[1]);
	if( !dataPoint::instance()->DPKin.isWithinDP() ) return 0;
	return std::abs(evaluate());
}
double AmpRelBreitWignerRes::integral() const{
	double norm=0;
	double res, err;
//	const gsl_rng_type *T;
//	gsl_rng *r;

	size_t d=2;//dimension of function
	//set limits
//	double xLimit_low[2] = {0.9,0.9};
//	double xLimit_high[2] = {2.,2.};
//	gsl_monte_function G = {&eval,2,0};
//	size_t calls = 500000;
//	gsl_rng_env_setup ();
//	T = gsl_rng_default;
//	r=gsl_rng_alloc(T);

//	gsl_monte_miser_state *s = gsl_monte_miser_alloc (d);
//	gsl_monte_miser_integrate (&G, , xLimit_low, xLimit_high, 2, calls, r,s,&res, &err);
//	gsl_monte_miser_free(s);
	return norm;
}
