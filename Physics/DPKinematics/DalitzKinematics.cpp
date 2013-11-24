//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff - initial API
//     Mathias Michel - kinematic functions
//-------------------------------------------------------------------------------

#include <stdlib.h>
#include "Physics/DPKinematics/DalitzKinematics.hpp"
//#include "Physics/DPKinematics/DataPoint.hpp"
#include "Core/PhysConst.hpp"

#include <stdlib.h>
#include "gsl/gsl_monte.h"
#include "gsl/gsl_monte_plain.h"
#include "gsl/gsl_monte_miser.h"
#include "gsl/gsl_monte_vegas.h"
DalitzKinematics* DalitzKinematics::inst = NULL;

double DalitzKinematics::calcHelicityAngle(double invMassSq23, double invMassSq13, double M, double m1, double m2, double m3){
	/*
	 * Calculated the helicity angle of inv. mass of particles 2 and 3.
	 * The angle is measured between particle 1 and 2! This is our definition!
	 * The cosine of the angle is returned.
	 *
	 */
	double invMassSq12= getThirdVariableSq(invMassSq23,invMassSq13);
	double s = invMassSq23;
	double t = invMassSq13;
	double u = invMassSq12;
	double qSq = (s-(M+m1)*(M+m1))*(s-(M-m1)*(M-m1))/(4*s);
	double qPrimeSq = (s-(m2+m3)*(m2+m3))*(s-(m2-m3)*(m2-m3))/(4*s);
	double cosAngle = ( s*(t-u)+(M*M-m1*m1)*(m2*m2-m3*m3) )/(4*s*sqrt(qSq)*sqrt(qPrimeSq));

	return cosAngle;
}
DalitzKinematics::DalitzKinematics(std::string _nameMother, std::string _name1, std::string _name2, std::string _name3):
																		Br(0.0), nameMother(_nameMother), name1(_name1), name2(_name2), name3(_name3)
{
	M = PhysConst::instance()->getMass(_nameMother);
	m1 = PhysConst::instance()->getMass(_name1);
	m2 = PhysConst::instance()->getMass(_name2);
	m3 = PhysConst::instance()->getMass(_name3);
	spinM = PhysConst::instance()->getJ(_nameMother);
	spin1 = PhysConst::instance()->getJ(_name1);
	spin2 = PhysConst::instance()->getJ(_name2);
	spin3 = PhysConst::instance()->getJ(_name3);

	if(M==-999 || m1==-999|| m2==-999|| m3==-999) {
		std::cout<<"DPKinematics: ERROR: masses not set! EXIT!"<<std::endl;
		exit(1);
	}
	init();
};
DalitzKinematics::DalitzKinematics(double _M, double _Br, double _m1, double _m2, double _m3,
		std::string _nameMother, std::string _name1, std::string _name2, std::string _name3):
				M(_M), Br(_Br), m1(_m1), m2(_m2), m3(_m3),
				nameMother(_nameMother), name1(_name1), name2(_name2), name3(_name3)
{

	spinM = PhysConst::instance()->getJ(_nameMother);
	spin1 = PhysConst::instance()->getJ(_name1);
	spin2 = PhysConst::instance()->getJ(_name2);
	spin3 = PhysConst::instance()->getJ(_name3);
	init();
};
DalitzKinematics::DalitzKinematics(const DalitzKinematics& other):
				M(other.M), spinM(other.spinM), Br(other.Br),
				m1(other.m1), m2(other.m2), m3(other.m3),
				spin1(other.spin1), spin2(other.spin2), spin3(other.spin3),
				name1(other.name1), name2(other.name2), name3(other.name3)
{
	init();
};
void DalitzKinematics::init(){
	m23_sq_min= ((m2+m3)*(m2+m3));
	m23_sq_max=((M-m1)*(M-m1));
	m13_sq_min=((m1+m3)*(m1+m3));
	m13_sq_max=((M-m2)*(M-m2));
	m12_sq_min=((m1+m2)*(m1+m2));
	m12_sq_max=((M-m3)*(M-m3));
	m23_min=((m2+m3)); m23_max=((M-m1));
	m13_min=((m1+m3)); m13_max=((M-m2));
	m12_min=((m1+m2)); m12_max=((M-m3));
	_DPareaCalculated=0;

};
//! returns 1 if point is within PHSP otherwise 0
double phspFunc(double* x, size_t dim, void* param) {
	if(dim!=2) return 0;
	DalitzKinematics* kin =	static_cast<DalitzKinematics*>(param);
	double m12sq = kin->getThirdVariableSq(x[0],x[1]);
	//use MC integration of a step function to obtain the DP area
	if( kin->isWithinDP(x[1],x[0],m12sq) ) return 1.0;
	return 0.0;
};

double DalitzKinematics::getDParea(){
	if(!_DPareaCalculated) calcDParea();
	return _DParea;
}
void DalitzKinematics::calcDParea(){

	//	std::cout<<"DPKinematics: DEBUG: calculating dalitz plot area"<<std::endl;
	size_t dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	double xLimit_low[2] = {m13_sq_min,m23_sq_min};
	double xLimit_high[2] = {m13_sq_max,m23_sq_max};

	size_t calls = 100000;
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator

	gsl_monte_function F = {&phspFunc,dim, const_cast<DalitzKinematics*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&F, xLimit_low, xLimit_high, 2, calls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	std::cout<<"DPKinematics: INFO: Area of dalitz plot form MC integration: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res<<std::endl;

	_DParea=res;
	_DPareaCalculated=1;
	return;
}
unsigned int DalitzKinematics::getSpin(unsigned int num){
	switch(num){
	case 0: return spinM;
	case 1: return spin1;
	case 2: return spin2;
	case 3: return spin3;
	}
	std::cout<<"DPKinematics: ERROR: wrong particle requested!"<<std::endl;
	exit(1);
	return -999;
}
unsigned int DalitzKinematics::getSpin(std::string name){
	if(name==nameMother) return spinM;
	if(name==name1) return spin1;
	if(name==name2) return spin2;
	if(name==name3) return spin3;
	std::cout<<"DPKinematics: ERROR: wrong particle requested!"<<std::endl;
	exit(1);
	return -999;
}

double DalitzKinematics::getMass(unsigned int num){
	switch(num){
	case 0: return M;
	case 1: return m1;
	case 2: return m2;
	case 3: return m3;
	}
	std::cout<<"DPKinematics: ERROR: wrong particle requested!"<<std::endl;
	exit(1);
	return -999;
}
double DalitzKinematics::getMass(std::string name){
	if(name==nameMother) return M;
	if(name==name1) return m1;
	if(name==name2) return m2;
	if(name==name3) return m3;
	std::cout<<"DPKinematics: ERROR: wrong particle requested!"<<std::endl;
	exit(1);
	return -999;
}

double DalitzKinematics::getThirdVariableSq(double invmass1sq, double invmass2sq) const{
	/*!
	 * calculates 3rd invariant mass from the other inv. masses.
	 */
	//	return sqrt(M*M+m1*m1+m2*m2+m3*m3-invmass1-invmass2);
	return (M*M+m1*m1+m2*m2+m3*m3-invmass1sq-invmass2sq);
}

//bool DPKinematics::isWithinDP() const{
//	/*!
//	 * checks if phase space point lies within the kinematically
//	 * allowed region. Point is taken from dataPoint singleton.
//	 */
//	static dataPoint* point = dataPoint::instance();
//	double s1 = point->getMsq(2,3);
//	double s2 = point->getMsq(1,3);
//	double s3 = point->getMsq(1,2);
//	return isWithinDP(s1,s2,s3);
//}

bool DalitzKinematics::isWithinDP(double m23sq, double m13sq, double m12sq) const{
	/*!
	 * \brief checks if phase space point lies within the kinematically allowed region.
	 * \param m23 invariant mass of particles 2 and 3
	 * \param m13 invariant mass of particles 2 and 3
	 * \param m12 invariant mass of particles 2 and 3
	 *
	 */
	//mostly copied from Laura++
	double e3Cms23 = (m23sq - m2*m2 + m3*m3)/(2.0*sqrt(m23sq)); //energy of part3 in 23 rest frame
	double p3Cms23 = sqrt(-(m3*m3-e3Cms23*e3Cms23)); //momentum of part3 in 23 rest frame
	double e1Cms23 = (M*M - m23sq - m1*m1)/(2.0*sqrt(m23sq)); //energy of part1 in 23 rest frame
	double p1Cms23 = sqrt(-(m1*m1-e1Cms23*e1Cms23)); //momentum of part1 in 23 rest frame

	double term = 2.0*e3Cms23*e1Cms23+ m1*m1 + m3*m3;

	double _m13SqLocMin = term - 2.0*p3Cms23*p1Cms23;
	double _m13SqLocMax = term + 2.0*p3Cms23*p1Cms23;

	if(m13sq >= _m13SqLocMin && m13sq <= _m13SqLocMax) return 1;
	return 0;
}

double DalitzKinematics::lambda(double x, double y, double z)const{
	return x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z;
}

double DalitzKinematics::s2min(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );

	double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) - lamterm )/(2.*s1);

	return result;
}

double DalitzKinematics::s2max(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );

	double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) + lamterm )/(2.*s1);

	return result;
}

double DalitzKinematics::s3min(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );

	double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m2*m2) - lamterm )/(2.*s1);

	return result;
}

double DalitzKinematics::s3max(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );

	double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m3*m3) + lamterm )/(2.*s1);

	return result;
}

double DalitzKinematics::s1min(double s2, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m3*m3, m1*m1) );

	double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) - lamterm )/(2.*s2);

	return result;
}

double DalitzKinematics::s1max(double s2, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m1*m1, m3*m3) );

	double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) + lamterm )/(2.*s2);

	return result;
}
