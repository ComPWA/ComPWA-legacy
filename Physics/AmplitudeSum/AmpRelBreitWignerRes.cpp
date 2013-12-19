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

//double integral(AmpRelBreitWignerRes* amp) {
//	double norm=0;
//	double res, err;
//	const gsl_rng_type *T;
//	gsl_rng *r;
//
//	size_t d=2;//dimension of function
//	//set limits
//	double xLimit_low[2] = {0.9,0.9};
//	double xLimit_high[2] = {2.,2.};
//	boost::function<unsigned int (AmpRelBreitWignerRes*,double*,size_t,void*)> sp;
//	sp = &AmpRelBreitWignerRes::eval;
//	gsl_monte_function G = {sp(amp),2,0};
//	size_t calls = 500000;
//	gsl_rng_env_setup ();
//	T = gsl_rng_default;
//	r=gsl_rng_alloc(T);
//
//	gsl_monte_miser_state *s = gsl_monte_miser_alloc (d);
//	//	gsl_monte_miser_integrate (&G, , xLimit_low, xLimit_high, 2, calls, r,s,&res, &err);
//	gsl_monte_miser_free(s);
//	return norm;
//}


AmpRelBreitWignerRes::AmpRelBreitWignerRes(const char *name,
		DoubleParameter& resMass, DoubleParameter& resWidth,
		double& mesonRadius, //  meson radius
		int subSys,
		int resSpin, int m, int n
) :
AmpAbsDynamicalFunction(name),
AmpKinematics(resMass, subSys, resSpin, m, n, AmpKinematics::barrierType(BWPrime), mesonRadius, 1.5),
_resWidth(resWidth),
//_wignerD(name, resSpin, m, n, subSys)
_wignerD(subSys, resSpin)
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
//	_wignerD.setDecayMasses(ma, mb, mc, M);
	return;
}
std::complex<double> AmpRelBreitWignerRes::evaluateAmp(dataPoint& point) const {
	if(_ma == -999 || _mb == -999 ||_mc == -999 ||_M == -999){
		std::cout<<"Masses of decay products not set!"<<std::endl;
		return 0;
	}
	std::complex<double> result;
//	double m = dataPoint::instance()->getM(_subSys);
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double m23sq = point.getVal("m23sq");
	double m13sq = point.getVal("m13sq");
	double m12sq = kin->getThirdVariableSq(m23sq,m13sq);
	double m = -999;
	switch(_subSys){
	case 3: m=sqrt(m12sq); break;
	case 4: m=sqrt(m13sq); break;
	case 5: m=sqrt(m23sq); break;
	}
	//	double spinTerm = evaluateWignerD(); //spinTerm =1;
	double BLWeiss2 = BLres2(m);
	double qTerm = pow(q(m) / q0(), 2.*_spin + 1.);
	double Gamma0 = _resWidth.GetValue();
	double GammaV = Gamma0 * qTerm * (_mR / m) * BLWeiss2;
	std::complex<double> denom(_mR*_mR - m*m, -_mR * GammaV);

	result = std::complex<double>( _norm ) / denom; //Laura++ (old) definition
//	result = std::complex<double>(_norm*sqrt(BLWeiss2)) / denom;
//	result = std::complex<double>( _norm*sqrt(BLWeiss2)*sqrt(BLmother2(m)) ) / denom; //Laura++ (new) definition

	if(result.real()!=result.real()) {std::cout << "RE part NAN" << std::endl;return 0;}
	if(result.imag()!=result.imag()) {std::cout << "IM part NAN" << std::endl; return 0;}
	return result;
}
std::complex<double> AmpRelBreitWignerRes::evaluate(dataPoint& point) const {
//	std::cout<<evaluateAmp()<<" "<<evaluateWignerD()<<std::endl;
	unsigned int twoJplusOne = (2*_spin+1);
//	return evaluateAmp()*evaluateWignerD(); //DEBUG
//	return evaluateAmp()*(double)twoJplusOne; //DEBUG
	return evaluateAmp(point)*evaluateWignerD(point)*(double)twoJplusOne;
}
/*std::complex<double> AmpRelBreitWignerRes::evaluateTree(const ParameterList& paras) const {
    double Gamma0, GammaV;
    double m0 = double(paras.GetParameterValue("m0_"+_name));
    double m  = double(paras.GetParameterValue("x_"+_name));
    double ma = double(paras.GetParameterValue("ma_"+_name));
    double mb = double(paras.GetParameterValue("mb_"+_name));
    unsigned int spin = double(paras.GetParameterValue("spin_"+_name));
    double d = double(paras.GetParameterValue("d_"+_name));
    double norm = double(paras.GetParameterValue("norm_"+_name));

    Gamma0 = double(paras.GetParameterValue("resWidth_"+_name));
   // GammaV = Gamma0 * (m0 / m) * pow(q(ma,mb,m) / q0(ma,mb,m0), 2.*spin + 1.)  * BLprime2(ma,mb,m0,m,d,spin);
	double BLWeiss2 = BLres2(m);
	double qTerm = pow(q(m) / q0(), 2.*spin + 1.);
	//double Gamma0 = _resWidth.GetValue();
	GammaV = Gamma0 * qTerm * (m0 / m) * BLWeiss2;
	std::complex<double> denom(m0*m0 - m*m, -m0 * GammaV);

    //RooComplex denom = RooComplex(m0*m0 - m*m, -m0 * GammaV);
    //RooComplex res = RooComplex(m0 * Gamma0) / denom; //TODO: use same functions as over evaluates

    std::complex<double> result = std::complex<double>( norm ) / denom;
    //std::shared_ptr<ComplexParameter> bw(new ComplexParameter("relBW of "+name, result));
    return result*_wignerD.evaluateTree(paras,_name);
}*/
