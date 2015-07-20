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

#include <stdlib.h>
#include "gsl/gsl_monte.h"
#include "gsl/gsl_monte_plain.h"
#include "gsl/gsl_monte_miser.h"
#include "gsl/gsl_monte_vegas.h"

#include "Core/PhysConst.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"

AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(const char *name,
		std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
		std::shared_ptr<DoubleParameter> mass, int subSys, Spin spin, Spin m, Spin n,
		std::shared_ptr<DoubleParameter> mesonR, //  meson radius
		std::shared_ptr<DoubleParameter> motherR, //  mother radius
		int nCalls, normStyle nS) :
		_name(name), _mag(mag), _phase(phase), _mass(mass), _subSys(subSys), _spin(spin),
		_m(m), _n(n), _mesonRadius(mesonR), _motherRadius(motherR),  _nCalls(nCalls),
		_normStyle(nS), _norm(1.0), modified(1), _wignerD(subSys, spin)
{
	initialize();
}
AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(const char *name,
		std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
		std::shared_ptr<DoubleParameter> mass, int subSys, Spin spin, Spin m, Spin n,
		int nCalls, normStyle nS) :
				_name(name), _mag(mag), _phase(phase), _mass(mass), _subSys(subSys), _spin(spin),
				_m(m), _n(n),
				_mesonRadius(std::make_shared<DoubleParameter>(name, 1.0)),
				_motherRadius(std::make_shared<DoubleParameter>(name, 1.0)),
				_nCalls(nCalls), _normStyle(nS), _norm(1.0), modified(1), _wignerD(subSys, spin)
{
	initialize();
}

void AmpAbsDynamicalFunction::initialize()
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

AmpAbsDynamicalFunction::~AmpAbsDynamicalFunction() 
{
}

std::complex<double> AmpAbsDynamicalFunction::evaluate(dataPoint& point){
	if(_mass->GetValue() != tmp_mass){
		SetModified();
		tmp_mass = _mass->GetValue();
	}

	double a = std::abs(_mag->GetValue());
	double phi = _phase->GetValue();
	std::complex<double> eiphi(a*cos(phi),a*sin(phi));
	std::complex<double> res = evaluateAmp(point);
	double ang = evaluateWignerD(point);
	double norm = GetNormalization();

	return (eiphi*norm*res*ang);
}

double evalAmp(double* x, size_t dim, void* param) {
	/* We need a wrapper here because a eval() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	dataPoint pp; pp.setVal(0,x[1]);pp.setVal(1,x[0]);
	if( !kin->isWithinPhsp(pp) ) return 0;//only integrate over phase space
	std::complex<double> res = static_cast<AmpAbsDynamicalFunction*>(param)->evaluateAmp(pp);
	return ( std::norm(res) ); //integrate over |F|^2
}

double AmpAbsDynamicalFunction::integral(){
	size_t dim=2;
	double res=0.0, err=0.0;

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function F = {&evalAmp,dim, const_cast<AmpAbsDynamicalFunction*> (this)};
	//	gsl_monte_function F = {&twoDimGaussian,dim, new int()};//using test function; result should be 1

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&F, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
//	BOOST_LOG_TRIVIAL(debug)<<"AmpAbsDynamicalFunction::integral() Integration result for |"
//			<<_name<<"|^2: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

double AmpAbsDynamicalFunction::GetNormalization(){
	if(_norm<0) return 1.0; //normalization is disabled
	if(!modified) return _norm;
	_norm = 1/sqrt(integral());
	modified=0;
	return _norm;
}

double eval(double* x, size_t dim, void* param) {
	/* We need a wrapper here because evaluate() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	dataPoint pp; pp.setVal(0,x[1]);pp.setVal(1,x[0]);
	if( !kin->isWithinPhsp(pp) ) return 0;//only integrate over phase space
	std::complex<double> res = static_cast<AmpAbsDynamicalFunction*>(param)->evaluateAmp(pp);
	double ang = static_cast<AmpAbsDynamicalFunction*>(param)->evaluateWignerD(pp);
	double norm = static_cast<AmpAbsDynamicalFunction*>(param)->GetNormalization();
	return ( std::norm(res*ang*norm) ); //integrate over |F|^2
}

double AmpAbsDynamicalFunction::totalIntegral() const{
	BOOST_LOG_TRIVIAL(debug)<<"AmpAbsDynamicalFunction::totalIntegral() calculating integral of "<<_name<<" !";
	size_t dim=2;
	double res=0.0, err=0.0;

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function F = {&eval,dim, const_cast<AmpAbsDynamicalFunction*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&F, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	BOOST_LOG_TRIVIAL(debug)<<"AmpAbsDynamicalFunction::totalIntegral() result for |"<<_name<<"|^2: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

std::complex<double> AmpAbsDynamicalFunction::widthToCoupling(double mSq, double mR, double width,
		double ma, double mb, double spin, double mesonRadius)
{
	double sqrtS = sqrt(mSq);
	//calculate gammaA(s_R)
	double ffR = Kinematics::FormFactor(mR,ma,mb,spin,mesonRadius);
	std::complex<double> qR = Kinematics::qValue(mR,ma,mb);
	//calculate phsp factor
	std::complex<double> rho = Kinematics::phspFactor(sqrtS,ma,mb);
	std::complex<double> denom = std::pow(qR,spin)*ffR*sqrt(rho);
	std::complex<double> result = std::complex<double>(sqrt(mR*width), 0) / denom;
	return result;
}

std::complex<double> AmpAbsDynamicalFunction::couplingToWidth(double mSq, double mR, double g,
		double ma, double mb, double spin, double mesonRadius)
{
	double sqrtM = sqrt(mSq);
	//calculate gammaA(s_R)
	double ffR = Kinematics::FormFactor(mR,ma,mb,spin,mesonRadius);
	std::complex<double> qR = std::pow(Kinematics::qValue(mR,ma,mb),spin);
	std::complex<double> gammaA = ffR*qR;
	//calculate phsp factor
	std::complex<double> rho = Kinematics::phspFactor(sqrtM,ma,mb);
	std::complex<double> result = std::norm(gammaA)*g*g*rho/ mR;
	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"AmpKinematics::couplingToWidth() | NaN! mSq="<<mSq
				<<" mR="<<mR<<" g="<<g<<" ma="<<ma<<" mb="<<mb<<std::endl;
		std::cout<<qR<<" "<<gammaA<<" "<<rho<<" "<<g<<std::endl;
	}
	return result;
}
double twoDimGaussian(double* z, size_t dim, void *param){
	if(dim!=2) return 0;
	/* test environment for numeric integration:
	 * 	Calculating integral of normalized gaussian:
	 * 	f(x,y) = A exp( - (x-x0)^2/(2*sigmaX^2) + (y-y0)^2/(2*sigmaY^2)
	 * 	with A=1/(2*pi*sigmaX*sigmaY) this function is normalized to 1
	 */
	double x = z[0]; double y = z[1];
	//mean and width need to be adjusted according to final state kinematics
	double x0=1.1, y0=1.1; //mean
	double sigmaX=0.01, sigmaY=0.01; //width
	double pi = PhysConst::instance()->getConstValue("Pi");

	double result = exp( -(x-x0)*(x-x0)/(2*sigmaX*sigmaX) - (y-y0)*(y-y0)/(2*sigmaY*sigmaY) );
	result/=2*pi*sigmaY*sigmaX;
	return result;
}
