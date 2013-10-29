/*
 * DPKinematics.cpp

 *
 *  Created on: Oct 16, 2013
 *      Author: weidenka
 */

#include "Physics/DPKinematics/DPKinematics.hpp"
#include "Physics/DPKinematics/DataPoint.hpp"
#include "Physics/DPKinematics/PhysConst.hpp"

void DPKinematics::init(){
	m23_sq_min= ((m2+m3)*(m2+m3));
	m23_sq_max=((M-m1)*(M-m1));
	m13_sq_min=((m1+m3)*(m1+m3));
	m13_sq_max=((M-m2)*(M-m2));
	m12_sq_min=((m1+m2)*(m1+m2));
	m12_sq_max=((M-m3)*(M-m3));
	m23_min=((m2+m3)); m23_max=((M-m1));
	m13_min=((m1+m3)); m13_max=((M-m2));
	m12_min=((m1+m2)); m12_max=((M-m3));
};
DPKinematics::DPKinematics(std::string _nameMother, std::string _name1, std::string _name2, std::string _name3):
						Br(0.0), nameMother(_nameMother), name1(_name1), name2(_name2), name3(_name3)
{
	M = PhysConst::instance()->getMass(_nameMother);
	m1 = PhysConst::instance()->getMass(_name1);
	m2 = PhysConst::instance()->getMass(_name2);
	m3 = PhysConst::instance()->getMass(_name3);
//	std::cout<<M<< " "<<m1<<" " <<m2<<" " <<m3<<std::endl;
	init();
};
DPKinematics::DPKinematics(double _M, double _Br, double _m1, double _m2, double _m3, std::string _name1, std::string _name2, std::string _name3):
						M(_M), Br(_Br), m1(_m1), m2(_m2), m3(_m3), name1(_name1), name2(_name2), name3(_name3)
{
	init();
};
DPKinematics::DPKinematics(const DPKinematics& other):
						M(other.M), Br(other.Br),
						m1(other.m1), m2(other.m2), m3(other.m3),
						name1(other.name1), name2(other.name2), name3(other.name3)
{
	init();
};
double DPKinematics::getThirdVariable(double invmass1, double invmass2) const{
	return sqrt(M*M+m1*m1+m2*m2+m3*m3-invmass1*invmass1-invmass2*invmass2);
}
bool DPKinematics::isWithinDP() const{
	/*!
	 * checks if phase space point lies within the kinematically
	 * allowed region. Point is taken from dataPoint singleton.
	 */
	static dataPoint* point = dataPoint::instance();
	double s1 = point->getMsq(2,3);
	double s2 = point->getMsq(1,3);
	double s3 = point->getMsq(1,2);
	return isWithinDP(s1,s2,s3);
}

bool DPKinematics::isWithinDP(double m23, double m13, double m12) const{
	/*!
	 * \brief checks if phase space point lies within the kinematically allowed region.
	 * \param m23 invariant mass of particles 2 and 3
	 * \param m13 invariant mass of particles 2 and 3
	 * \param m12 invariant mass of particles 2 and 3
	 *
	 */
	//mostly copied from Laura++
	double e3Cms23 = (m23*m23 - m2*m2 + m3*m3)/(2.0*m23); //energy of part3 in 23 rest frame
	double p3Cms23 = sqrt(-(m3*m3-e3Cms23*e3Cms23)); //momentum of part3 in 23 rest frame
	double e1Cms23 = (M*M - m23*m23 - m1*m1)/(2.0*m23); //energy of part1 in 23 rest frame
	double p1Cms23 = sqrt(-(m1*m1-e1Cms23*e1Cms23)); //momentum of part1 in 23 rest frame

	double term = 2.0*e3Cms23*e1Cms23+ m1*m1 + m3*m3;

	double _m13SqLocMin = term - 2.0*p3Cms23*p1Cms23;
	double _m13SqLocMax = term + 2.0*p3Cms23*p1Cms23;

	if(m13*m13 >= _m13SqLocMin && m13*m13 <= _m13SqLocMax) return 1;
	return 0;
}

double DPKinematics::lambda(double x, double y, double z)const{
	return x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z;
}

double DPKinematics::s2min(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );

	double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) - lamterm )/(2.*s1);

	return result;
}

double DPKinematics::s2max(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );

	double result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) + lamterm )/(2.*s1);

	return result;
}

double DPKinematics::s3min(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );

	double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m2*m2) - lamterm )/(2.*s1);

	return result;
}

double DPKinematics::s3max(double s1, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );

	double result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m3*m3) + lamterm )/(2.*s1);

	return result;
}

double DPKinematics::s1min(double s2, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m3*m3, m1*m1) );

	double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) - lamterm )/(2.*s2);

	return result;
}

double DPKinematics::s1max(double s2, double m0, double m1, double m2, double m3)const
{
	double s      = m0*m0;
	double lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m1*m1, m3*m3) );

	double result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) + lamterm )/(2.*s2);

	return result;
}
