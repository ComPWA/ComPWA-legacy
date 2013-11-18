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

//! DPKinematics stores kinematic information of a dalitz plot
/*!
 * @file DPKinematics.hpp
 *\class DPKinematics
 * 			DPKinematics stores particle masses of a process and provides some functions which calculate kinematic
 * 			quantities.
 */
#ifndef DPKINEMATICS_HPP_
#define DPKINEMATICS_HPP_

#include <iostream>
class DPpoint;
class dataPoint;

class DPKinematics
{
private:
	bool _DPareaCalculated;
	double _DParea;
	//! calculated dalitz plot area for the given kinematics
	void calcDParea();
//	std::shared_ptr<DalitzEfficiency> eff;
public:
	DPKinematics(){};
	void init();
	//! constructor access particles by name
	DPKinematics(std::string _nameMother, std::string _name1, std::string _name2, std::string _name3);
	//! constructor with particle masses and names
	DPKinematics(double _M, double _Br, double _m1, double _m2, double _m3,
			std::string _nameMother, std::string _name1, std::string _name2, std::string _name3);
	//! Copy constructor
	DPKinematics(const DPKinematics& other);
	//! calculated the helicity angle
	double calcHelicityAngle(double invM1sq, double invM2sq, double M, double ma, double mb, double mc);
	//! Calculates third dalitz plot variable, e.g f(s1,s2)=s3
	double getThirdVariableSq(double, double) const;
	//! checks of data point is within phase space boundaries, data point provided by dataPoint
	bool isWithinDP() const;
	//! checks of data point is within phase space boundaries
	bool isWithinDP(double m23, double m13, double m12=0) const;
	//! returns the dalitz plot area for the given kinematics
	double getDParea();

	double lambda(double x, double y, double z)const;
	double s2min(double s1, double m0, double m1, double m2, double m3)const;
	double s2max(double s1, double m0, double m1, double m2, double m3)const;
	double s3min(double s1, double m0, double m1, double m2, double m3)const;
	double s3max(double s1, double m0, double m1, double m2, double m3)const;
	double s1min(double s2, double m0, double m1, double m2, double m3)const;
	double s1max(double s2, double m0, double m1, double m2, double m3)const;
	double s2min(double s1)const { return s2min(s1,M,m1,m2,m3); };
	double s2max(double s1)const { return s2max(s1,M,m1,m2,m3); };
	double s3min(double s1)const { return s3min(s1,M,m1,m2,m3); };
	double s3max(double s1)const { return s3max(s1,M,m1,m2,m3); };
	double s1min(double s2)const { return s1min(s2,M,m1,m2,m3); };
	double s1max(double s2)const { return s1max(s2,M,m1,m2,m3); };

	//! get mass of paticles
	double getMass(unsigned int num);
	//! get mass of paticles
	double getMass(std::string name);
	//! get spin of decaying particle
	unsigned int getSpin(unsigned int num);
	//! get spin of particles
	unsigned int getSpin(std::string name);

	//! mass of decaying particle
	double M; unsigned int spinM;
	//! width of decaying particle
	double Br;
	//! masses of final state particles
	double m1; unsigned int spin1;
	double m2; unsigned int spin2;
	double m3; unsigned int spin3;

	//! names of particles
	std::string nameMother, name1, name2, name3;

	double m23_sq_min, m23_sq_max;
	double m13_sq_min, m13_sq_max;
	double m12_sq_min, m12_sq_max;
	double m23_min, m23_max;
	double m13_min, m13_max;
	double m12_min, m12_max;

};

#endif
