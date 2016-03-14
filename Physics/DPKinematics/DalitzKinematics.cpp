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
#include "Core/DataPoint.hpp"

#include <stdlib.h>
#include "gsl/gsl_monte.h"
#include "gsl/gsl_monte_plain.h"
#include "gsl/gsl_monte_miser.h"
#include "gsl/gsl_monte_vegas.h"
#include "gsl/gsl_sf_legendre.h"

namespace ComPWA {
namespace Physics {
namespace DPKinematics {

DalitzKinematics::DalitzKinematics(std::string _nameMother,
		std::string _name1, std::string _name2, std::string _name3) :
		Br(0.0), nameMother(_nameMother),
		name1(_name1), name2(_name2), name3(_name3), massIdsSet(false)
{
	number_of_particles_ = 3;
	M = ComPWA::PhysConst::Instance().findParticle(_nameMother).mass_;
	m1 = ComPWA::PhysConst::Instance().findParticle(_name1).mass_;
	m2 = ComPWA::PhysConst::Instance().findParticle(_name2).mass_;
	m3 = ComPWA::PhysConst::Instance().findParticle(_name3).mass_;
	spinM = ComPWA::PhysConst::Instance().findParticle(_nameMother).spin_.J_numerator_;
	spin1 = ComPWA::PhysConst::Instance().findParticle(_name1).spin_.J_numerator_;
	spin2 = ComPWA::PhysConst::Instance().findParticle(_name2).spin_.J_numerator_;
	spin3 = ComPWA::PhysConst::Instance().findParticle(_name3).spin_.J_numerator_;

	BOOST_LOG_TRIVIAL(info) << " DalitzKinematics::DalitzKinematics() | Setting up decay "
			<<_nameMother<<"->"<<_name1<<" "<<_name2<<" "<<_name3;
	init();
};
DalitzKinematics::DalitzKinematics(double _M, double _Br, double _m1, double _m2, double _m3,
		std::string _nameMother, std::string _name1, std::string _name2, std::string _name3):
		M(_M), Br(_Br), m1(_m1), m2(_m2), m3(_m3),
		nameMother(_nameMother), name1(_name1), name2(_name2), name3(_name3), massIdsSet(false)
{

	number_of_particles_ = 3;
  spinM = ComPWA::PhysConst::Instance().findParticle(_nameMother).spin_.J_numerator_;
  spin1 = ComPWA::PhysConst::Instance().findParticle(_name1).spin_.J_numerator_;
  spin2 = ComPWA::PhysConst::Instance().findParticle(_name2).spin_.J_numerator_;
  spin3 = ComPWA::PhysConst::Instance().findParticle(_name3).spin_.J_numerator_;
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
	Msq = M*M;
	mSq1 = m1*m1;
	mSq2 = m2*m2;
	mSq3 = m3*m3;
	mSq4 = m4*m4;
	
	variable_names_.push_back("m23sq");
	variable_names_.push_back("m13sq");
	//	varNames.push_back("m12sq");
	
	BOOST_LOG_TRIVIAL(debug) << "PHSP boundaries:";
	BOOST_LOG_TRIVIAL(debug) << " "<<m23_min<<" < m23 < "<<m23_max<< " | "<<m23_sq_min<<" < m23sq < "<<m23_sq_max;
	BOOST_LOG_TRIVIAL(debug) << " "<<m13_min<<" < m13 < "<<m13_max<< " | "<<m13_sq_min<<" < m13sq < "<<m13_sq_max;
	BOOST_LOG_TRIVIAL(debug) << " "<<m12_min<<" < m12 < "<<m12_max<< " | "<<m12_sq_min<<" < m12sq < "<<m12_sq_max;

	return;
};
double DalitzKinematics::calculateMoments(unsigned int sys, dataPoint& point, unsigned int n, unsigned int m){
	double angle = helicityAngle(sys,point);
	if(angle < -1 || angle > 1 ) {
		BOOST_LOG_TRIVIAL(error) << "DalitzKinematics::calculateMoments() angle out of range! "<<angle;
		return -999;
	}
	//	double val = gsl_sf_legendre_Plm(n,m,angle);
	double val = gsl_sf_legendre_sphPlm(n,m,angle);//normalized! - which one is correct?
	//	std::cout<<angle<< " "<<val<<std::endl;
	return val;
}
double DalitzKinematics::getMin(std::string name){
	if(name=="m13sq") return m13_sq_min;
	if(name=="m23sq") return m23_sq_min;
	if(name=="m12sq") return m12_sq_min;
	return -1;
}
double DalitzKinematics::getMax(std::string name){
	if(name=="m13sq") return m13_sq_max;
	if(name=="m23sq") return m23_sq_max;
	if(name=="m12sq") return m12_sq_max;
	return -1;
}

void DalitzKinematics::translateEventToDataPoint(const Event& ev, dataPoint& point) const {
	Particle part1 = ev.getParticle(0);
	Particle part2 = ev.getParticle(1);
	Particle part3 = ev.getParticle(2);
	double m23sq = Particle::invariantMass(part2,part3);
	double m13sq = Particle::invariantMass(part1,part3);
	//	point.setVal("m13sq",m13sq); point.setVal("m23sq",m23sq);
	point.setVal(1,m13sq); point.setVal(0,m23sq);
	return;
}
double DalitzKinematics::invMassMax(unsigned int sys, unsigned int sys2, double invMass_sys) const {
	unsigned int part1, part2;
	double Mfirst, Msecond;
	switch(sys){
	case 5: part1=2; part2=3; Mfirst=m2; Msecond=m3; break;
	case 4: part1=3; part2=1;Mfirst=m3; Msecond=m1; break;
	case 3: part1=1; part2=2;Mfirst=m1; Msecond=m2; break;
	}
	double Efirst=eiCms(part1,sys2,invMass_sys);
	double Esecond=eiCms(part2,sys2,invMass_sys);

	double max=(Efirst+Esecond)*(Efirst+Esecond)-(sqrt(Efirst*Efirst-Mfirst*Mfirst)-sqrt(Esecond*Esecond-Msecond*Msecond))*(sqrt(Efirst*Efirst-Mfirst*Mfirst)-sqrt(Esecond*Esecond-Msecond*Msecond));
	return max;
}
double DalitzKinematics::invMassMin(unsigned int sys, unsigned int sys2, double invMass_sys) const {
	unsigned int part1, part2;
	double Mfirst, Msecond;
	switch(sys){
	case 5: part1=2; part2=3; Mfirst=m2; Msecond=m3; break;
	case 4: part1=3; part2=1;Mfirst=m3; Msecond=m1; break;
	case 3: part1=1; part2=2;Mfirst=m1; Msecond=m2; break;
	}
	double Efirst=eiCms(part1,sys2,invMass_sys);
	double Esecond=eiCms(part2,sys2,invMass_sys);

	double min=(Efirst+Esecond)*(Efirst+Esecond)-(sqrt(Efirst*Efirst-Mfirst*Mfirst)+sqrt(Esecond*Esecond-Msecond*Msecond))*(sqrt(Efirst*Efirst-Mfirst*Mfirst)+sqrt(Esecond*Esecond-Msecond*Msecond));
	return min;
}
double DalitzKinematics::eiCms(unsigned int partId, unsigned int sys, double invMass_sys) const {
	double E;
	double m = sqrt(invMass_sys);
	switch(sys){
	case 3:
		switch(partId){
		case 1: E = (invMass_sys-mSq2+mSq1)/(2*m);break;
		case 2: E = (invMass_sys-mSq1+mSq2)/(2*m);break;
		case 3: E =-(invMass_sys-Msq+mSq3)/(2*m);	break;
		}
		break;
	case 4:
		switch(partId){
		case 1: E = (invMass_sys-mSq3+mSq1)/(2*m);break;
		case 2: E =-(invMass_sys-Msq+mSq2)/(2*m);break;
		case 3: E = (invMass_sys-mSq1+mSq3)/(2*m); break;
		}
		break;
	case 5:
		switch(partId){
		case 1: E =-(invMass_sys-Msq+mSq1)/(2*m);break;
		case 2: E = (invMass_sys-mSq3+mSq2)/(2*m);break;
		case 3: E = (invMass_sys-mSq2+mSq3)/(2*m); break;
		}
		break;
	}

	return E;
}
double DalitzKinematics::mimin(unsigned int i) const
{
	double result = 0.;

	if(i!=5&&i!=4&&i!=3) {
		BOOST_LOG_TRIVIAL(error)<<"DalitzKinematics: ERROR: Wrong subsystem! ";
		return -999;
	}
	switch (i)
	{
	case 5:  result = (m2+m3)*(m2+m3);
	break;
	case 4:  result = (m1+m3)*(m1+m3);
	break;
	case 3:  result = (m1+m2)*(m1+m2);
	break;
	}
	return result;
}

double DalitzKinematics::mimax(unsigned int i) const
{
	double result = 0.;

	if(i!=5&&i!=4&&i!=3) {
		BOOST_LOG_TRIVIAL(error)<<"DalitzKinematics: ERROR: Wrong subsystem! ";
		return -999;
	}
	switch (i)
	{
	case 5:  result = (M-m1)*(M-m1);
	break;
	case 4:  result = (M-m2)*(M-m2);
	break;
	case 3:  result = (M-m3)*(M-m3);
	break;
	}
	return result;
}

void DalitzKinematics::phspContour(unsigned int xsys,unsigned int ysys,unsigned int n, double* xcoord, double* ycoord)
{
	unsigned int num=n;
	if(num%2!=0) {
		num-=1;
		BOOST_LOG_TRIVIAL(info)<<"DalitzKinematics::phspContour(): setting size to a even number. Assure that the size of your arrays is "<<num*2+1<<"!";
	}
	//    g->Set(n*2+1);
	double massminl = mimin(xsys);
	double massmaxl = mimax(xsys);

	double binw = (massmaxl-massminl)/(double)(num);

	//	std::cout<<"==== x="<<xsys<<" y="<<ysys<<" "<<massminl<<" "<<massmaxl<<" "<<binw<<std::endl;
	double ymin,ymax,x;
	unsigned int i=0;

	for (; i<num; i++)
	{
		x = i*binw + massminl;
		while(x<=massminl) { x+=binw/100; }
		while(x>=massmaxl) { x-=binw/100; }
		ymin = invMassMin(ysys,xsys,x);
		ymax = invMassMax(ysys,xsys,x);

		xcoord[i]=x; ycoord[i]=ymin;
		xcoord[num*2-i]=x; ycoord[num*2-i]=ymax;
		//		std::cout<<x<< " "<<i<< " "<<ymin<<" "<<" "<<num*2-i<<" "<<ymax<<std::endl;
	}
	//adding last datapoint
	x = i*binw + massminl;
	while(x<=massminl) { x+=binw/100; }
	while(x>=massmaxl) { x-=binw/100; }

	ymin = invMassMin(ysys,xsys,x);
	ymax = invMassMax(ysys,xsys,x);

	//	std::cout<<x<< " "<<i<< " "<<ymin<<" "<<" "<<num*2-i<<" "<<ymax<<std::endl;
	xcoord[i]=x; ycoord[i]=ymin;
	return;
}

double DalitzKinematics::helicityAngle(unsigned int sys, const dataPoint& point){
	return helicityAngle(sys,point.getVal(0),point.getVal(1));
}
//double DalitzKinematics::calcHelicityAngle(unsigned int sys, dataPoint& point){
//	double cosTheta;
//	double cosTheta_test;
//	double m13sq = point.getVal(1);
//	double m23sq = point.getVal(0);
//	double m12sq = getThirdVariableSq(m23sq,m13sq);
//	//I have no idea how we should define the angles!
//	switch(sys){
//	case 3://angle versus K-
//		cosTheta = scatteringAngle(m12sq,m23sq,M,m3,m1,m2); break;//angle versus particle 2
//	case 4://angle versus K_S0
//		cosTheta = scatteringAngle(m13sq,m12sq,M,m2,m3,m1); break;//angle versus particle 1
//	case 5:
//		cosTheta = scatteringAngle(m23sq,m13sq,M,m1,m2,m3); break;//angle versus particle 3
//	default://angle versus K+
//		BOOST_LOG_TRIVIAL(error) <<"DalitzKinematics::calcHelicityAngle() : wrong subsystem selected!";
//		return 0.0;
//	}
//	//	    if(cosTheta>1.) cosTheta=1.;
//	//	    if(cosTheta<-1.) cosTheta=-1.;
//	return cosTheta;
//}

double DalitzKinematics::scatteringAngle(double s, double t, \
		double M, double mSpec, double mSecond, double m){
	double u = getThirdVariableSq(s,t);
	double qSq = (s-(M+mSpec)*(M+mSpec))*(s-(M-mSpec)*(M-mSpec))/(4*s);
	double qPrimeSq = (s-(mSecond+m)*(mSecond+m))*(s-(mSecond-m)*(mSecond-m))/(4*s);
	double cosAngle = ( s*(t-u)+(M*M-mSpec*mSpec)*(mSecond*mSecond-m*m) )/(4*s*sqrt(qSq)*sqrt(qPrimeSq));

	return cosAngle;
}
double DalitzKinematics::helicityAngle(unsigned int sys, double invMassSq23, double invMassSq13){
	double invMassSq12= getThirdVariableSq(invMassSq23,invMassSq13);
	double invMsqSys, invMsqSecond;
	double m, mSpec, eCms, eSpecCms, pCms, pSpecCms, cosAngle;
	switch(sys){
	case 3://angle versus particle 2
		m = m2; mSpec = m3; invMsqSys = invMassSq12; invMsqSecond=invMassSq23;
		eCms = eiCms(2,sys,invMsqSys);
		eSpecCms = eiCms(3,sys,invMsqSys);
		pCms = sqrt(eCms*eCms-m*m);
		pSpecCms = sqrt(eSpecCms*eSpecCms-mSpec*mSpec);
		break;
	case 4://angle versus particle 1
		m = m1; mSpec = m2; invMsqSys = invMassSq13; invMsqSecond=invMassSq12;
		eCms = eiCms(1,sys,invMsqSys);
		eSpecCms = eiCms(2,sys,invMsqSys);
		pCms = sqrt(eCms*eCms-m*m);
		pSpecCms = sqrt(eSpecCms*eSpecCms-mSpec*mSpec);
		break;
	case 5:
		//angle versus particle 3
//		m = m3; mSpec = m1; invMsqSys = invMassSq23; invMsqSecond=invMassSq13;
//		eCms = eiCms(3,sys,invMsqSys);
		//angle versus particle 2
		m = m2; mSpec = m1; invMsqSys = invMassSq23; invMsqSecond=invMassSq12;
		eCms = eiCms(2,sys,invMsqSys);
		eSpecCms = eiCms(1,sys,invMsqSys);
		pCms = sqrt(eCms*eCms-m*m);
		pSpecCms = sqrt(eSpecCms*eSpecCms-mSpec*mSpec);
		break;
	}
	cosAngle = -(invMsqSecond-m*m-mSpec*mSpec-2*eCms*eSpecCms)/(2*pCms*pSpecCms);

	return cosAngle;
}

//! returns 1 if point is within PHSP otherwise 0
double phspFunc(double* x, size_t dim, void* param) {
	if(dim!=2) return 0;
	DalitzKinematics* kin =	static_cast<DalitzKinematics*>(param);
	//	double m12sq = kin->getThirdVariableSq(x[0],x[1]);
	//use MC integration of a step function to obtain the DP area
	dataPoint pp; pp.setVal(0,x[1]);pp.setVal(1,x[0]);
	if( kin->isWithinPhsp(pp) ) return 1.0;
	return 0.0;
};

double DalitzKinematics::calculatePSArea() {
	size_t dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	double xLimit_low[2] = {m13_sq_min,m23_sq_min};
	double xLimit_high[2] = {m13_sq_max,m23_sq_max};

	size_t calls = 2000000;
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator

	gsl_monte_function F = {&phspFunc,dim, const_cast<DalitzKinematics*> (this)};

	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&F, xLimit_low, xLimit_high, 2, calls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	BOOST_LOG_TRIVIAL(debug)<<"DPKinematics::calcDParea() Dalitz plot area (MC integration): "
			<<"("<<res<<"+-"<<err<<") GeV^4 relAcc [%]: "<<100*err/res;

	return res;
}
unsigned int DalitzKinematics::getSpin(unsigned int num){
	switch(num){
	case 0: return spinM;
	case 1: return spin1;
	case 2: return spin2;
	case 3: return spin3;
	}
	throw std::runtime_error("DPKinematics::getSpin(int) | Wrong particle requested!");
	return -999;
}
unsigned int DalitzKinematics::getSpin(std::string name){
	if(name==nameMother) return spinM;
	if(name==name1) return spin1;
	if(name==name2) return spin2;
	if(name==name3) return spin3;
	throw std::runtime_error("DPKinematics::getSpin(string) | Wrong particle "+name+" requested!");
}

double DalitzKinematics::getMass(unsigned int num) const {
	switch(num){
	case 0: return M;
	case 1: return m1;
	case 2: return m2;
	case 3: return m3;
	}
	throw std::runtime_error("DPKinematics::getMass(int) | Wrong particle requested!");
}
double DalitzKinematics::getMass(std::string name) const {
	if(name==nameMother) return M;
	if(name==name1) return m1;
	if(name==name2) return m2;
	if(name==name3) return m3;
	throw std::runtime_error("DPKinematics::getMass(string) | Wrong particle "+name+" requested!");
}

double DalitzKinematics::getThirdVariableSq(double invmass1sq, double invmass2sq) const{
	/*!
	 * calculates 3rd invariant mass from the other inv. masses.
	 */
	//	return sqrt(M*M+m1*m1+m2*m2+m3*m3-invmass1-invmass2);
	return (Msq+mSq1+mSq2+mSq3-invmass1sq-invmass2sq);
}

bool DalitzKinematics::isWithinPhsp(const dataPoint& point) {
	if(!massIdsSet){
		id23 = point.getID("m23sq");
		id13 = point.getID("m13sq");
		massIdsSet = true;
	}
	double m23sq = point.getVal(id23);
	double m13sq = point.getVal(id13);
	double m12sq=getThirdVariableSq(m23sq,m13sq);

	bool c0=0; bool c1=0; bool c2=0;
	if(m13sq >= invMassMin(4,5,m23sq) && m13sq <= invMassMax(4,5,m23sq)) c0=1;
	if(m12sq >= invMassMin(3,5,m23sq) && m12sq <= invMassMax(3,5,m23sq)) c1=1;
	if(m23sq >= invMassMin(5,4,m13sq) && m23sq <= invMassMax(5,4,m13sq)) c2=1;
	if(c0 && c1 && c2) return 1;
	//if(c0 && c1) return 1;
	//if(c0) return 1;
	return 0;
}

//double DalitzKinematics::lambda(double x, double y, double z)const{
//	return x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z;
//}

//double DalitzKinematics::s2min(double s1, double m0, double m1, double m2, double m3) const {
//	double s      = m0*m0;
//	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );
//	return m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) - lamterm )/(2.*s1);
//}
//
//double DalitzKinematics::s2max(double s1, double m0, double m1, double m2, double m3)const
//{
//	double s      = m0*m0;
//	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );
//
//	double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) + lamterm )/(2.*s1);
//
//	return result;
//}
//
//double DalitzKinematics::s3min(double s1, double m0, double m1, double m2, double m3)const
//{
//	double s      = m0*m0;
//	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );
//
//	double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m2*m2) - lamterm )/(2.*s1);
//
//	return result;
//}
//
//double DalitzKinematics::s3max(double s1, double m0, double m1, double m2, double m3)const
//{
//	double s      = m0*m0;
//	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );
//
//	double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m3*m3) + lamterm )/(2.*s1);
//
//	return result;
//}
//
//double DalitzKinematics::s1min(double s2, double m0, double m1, double m2, double m3)const
//{
//	double s      = m0*m0;
//	double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m3*m3, m1*m1) );
//
//	double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) - lamterm )/(2.*s2);
//
//	return result;
//}
//
//double DalitzKinematics::s1max(double s2, double m0, double m1, double m2, double m3)const
//{
//	double s      = m0*m0;
//	double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m1*m1, m3*m3) );
//
//	double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) + lamterm )/(2.*s2);
//
//	return result;
//}

} /* namespace DPKinematics */
} /* namespace Physics */
} /* namespace ComPWA */
