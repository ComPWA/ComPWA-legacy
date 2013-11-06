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

#include "Physics/DPKinematics/PhysConst.hpp"
#include "Physics/DPKinematics/dataPoint.hpp"
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

double AmpAbsDynamicalFunction::absEvaluate(double x[], size_t dim) const {
	if(dim!=2) return 0;
	//set data point: we assume that x[0]=m13 and x[1]=m23
	double m12 = dataPoint::instance()->DPKin.getThirdVariable(sqrt(x[0]),sqrt(x[1]));
	dataPoint::instance()->setMsq(4,x[0]);
	dataPoint::instance()->setMsq(5,x[1]);
	dataPoint::instance()->setMsq(3,m12*m12);
	if( !dataPoint::instance()->DPKin.isWithinDP() ) return 0;//only integrate over phase space
	return std::abs(evaluate());
}
double evalWrapper(double* x, size_t dim, void* param) {
	/* We need a wrapper here because a eval() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	return static_cast<AmpAbsDynamicalFunction*>(param)->absEvaluate(x,dim);
};
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

double AmpAbsDynamicalFunction::integral() const{

//	std::cout<<"AmpRelBreitWignerRes: DEBUG: calculating integral of "<<_name<<" !"<<std::endl;
	size_t dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13 and x[1]=m23
	double xLimit_low[2] = {dataPoint::instance()->DPKin.m13_min,dataPoint::instance()->DPKin.m23_min};
	double xLimit_high[2] = {dataPoint::instance()->DPKin.m13_max,dataPoint::instance()->DPKin.m23_max};

//	double x[2]; x[0]=1.3; x[1]=1.301;//test values
//	double (*f)(double* x, size_t dim, void* param) = &evalWrapper;
//	std::cout<<"===== "<<f(x,2,const_cast<AmpRelBreitWignerRes*> (this))<<std::endl;

//	std::cout<<"AmpAbsDynamicalFunction: DEBUG: calculating integral in limits "\
//			<<xLimit_low[0]<<","<<xLimit_high[0]<<" and "<<xLimit_low[1]<<","<<xLimit_high[1]<<std::endl;

	size_t calls = 100000;
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator

	gsl_monte_function F = {&evalWrapper,dim, const_cast<AmpAbsDynamicalFunction*> (this)};
//	gsl_monte_function F = {&twoDimGaussian,dim, new int()};//using test function; result should be 1

//	gsl_monte_plain_state *s = gsl_monte_plain_alloc (dim);
//	gsl_monte_plain_integrate (&F, xLimit_low, xLimit_high, dim, calls, r,s,&res, &err);
//	gsl_monte_plain_free(s);
//	gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim);
//	gsl_monte_miser_integrate (&F, xLimit_low, xLimit_high, dim, calls, r,s,&res, &err);
//	gsl_monte_miser_free(s);

	/*	Choosing vegas algorithm here, because it is the most accurate:
	* 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	* 		 this should be sufficiency for most applications
	*/
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&F, xLimit_low, xLimit_high, 2, calls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	std::cout<<"AmpAbsDynamicalFunction: INFO: Integration result for "<<_name<<": "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res<<std::endl;

	return res;
}
