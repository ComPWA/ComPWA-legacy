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

AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(const char *name) : _name(name), _norm(1.0)
{

}

AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(const AmpAbsDynamicalFunction& other, const char* newname)
{
}

AmpAbsDynamicalFunction::~AmpAbsDynamicalFunction() 
{
}

//double AmpAbsDynamicalFunction::evaluate(double x[], size_t dim) {
//	if(dim!=2) return 0;
//	//set data point: we assume that x[0]=m13 and x[1]=m23
//	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
//	double m12sq = kin->getThirdVariableSq(x[0],x[1]);
//	dataPoint pp; pp.setVal("m23sq",x[1]);pp.setVal("m13sq",x[0]);
//	if( !kin->isWithinPhsp(pp) ) return 0;//only integrate over phase space
//	std::complex<double> res = evaluate(pp);
//	return ( std::abs(res)*std::abs(res) ); //integrate over |F|^2
//}
double evalWrapper(double* x, size_t dim, void* param) {
	/* We need a wrapper here because a eval() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
//	double m12sq = kin->getThirdVariableSq(x[0],x[1]);
	dataPoint pp; pp.setVal(0,x[1]);pp.setVal(1,x[0]);
	if( !kin->isWithinPhsp(pp) ) return 0;//only integrate over phase space
	std::complex<double> res = static_cast<AmpAbsDynamicalFunction*>(param)->evaluateAmp(pp);
	res = res * static_cast<AmpAbsDynamicalFunction*>(param)->evaluateWignerD(pp);
	return ( std::abs(res)*std::abs(res) ); //integrate over |F|^2
//	return static_cast<AmpAbsDynamicalFunction*>(param)->evaluate(x,dim);
};

double AmpAbsDynamicalFunction::integral() const{
	BOOST_LOG_TRIVIAL(debug)<<"AmpAbsDynamicalFunction::integral() calculating integral of "<<_name<<" !";
	size_t dim=2;
	double res=0.0, err=0.0;

		DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	size_t calls = 100000;
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function F = {&evalWrapper,dim, const_cast<AmpAbsDynamicalFunction*> (this)};
//	gsl_monte_function F = {&twoDimGaussian,dim, new int()};//using test function; result should be 1

	/*	Choosing vegas algorithm here, because it is the most accurate:
	* 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	* 		 this should be sufficiency for most applications
	*/
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&F, xLimit_low, xLimit_high, 2, calls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	BOOST_LOG_TRIVIAL(debug)<<"AmpAbsDynamicalFunction::integral() Integration result for |"<<_name<<"|^2: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
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
