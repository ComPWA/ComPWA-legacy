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
//		Klaus Goetzen - most of the kin. calculations
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
#include <vector>

#include <boost/log/trivial.hpp>
using namespace boost::log;

#include "Core/Kinematics.hpp"

class DalitzKinematics : public Kinematics
{
private:
//	static DalitzKinematics* inst;
	bool _DPareaCalculated;
	double _DParea;
	//! calculated dalitz plot area for the given kinematics
	void calcDParea();
	void init();

	//! Copy constructor
	DalitzKinematics(const DalitzKinematics& other);

	DalitzKinematics(){};
	//! constructor access particles by name
	DalitzKinematics(std::string _nameMother, std::string _name1, std::string _name2, std::string _name3);
	//! constructor with particle masses and names
	DalitzKinematics(double _M, double _Br, double _m1, double _m2, double _m3,
			std::string _nameMother, std::string _name1, std::string _name2, std::string _name3);

//	std::vector<std::string> varNames;


public:
	static Kinematics* createInstance(std::string _nameMother, std::string _name1, std::string _name2, std::string _name3){
		_inst = new DalitzKinematics(_nameMother, _name1, _name2,_name3);
		return _inst;
	}
//	static DalitzKinematics* instance(){
//		if(!inst) {
//			BOOST_LOG_TRIVIAL(fatal)<<"DPKinematics: ERROR instance not created first!";
//			return 0;
//		}
//		return inst;
//	};

	unsigned int sizeOfPhsp(){ return 3; }
//	std::vector<std::string> getVarNames(){return varNames;}
	void phspContour(unsigned int xsys,unsigned int ysys, unsigned int n, double* xcoord, double* ycoord);
	//! calculated the helicity angle
	double calcHelicityAngle(double invM1sq, double invM2sq, double M, double ma, double mb, double mc);
	//! Calculates third dalitz plot variable, e.g f(s1,s2)=s3
	double getThirdVariableSq(double, double) const;
	//! checks of data point is within phase space boundaries
	bool isWithinPhsp(const dataPoint &point) const;
	//! returns the dalitz plot area for the given kinematics
	double getDParea();

	//!maximum value for variable m23=5, m13=4, m12=3
	double mimax(unsigned int i) const;
	//!minimum value for variable m23=5, m13=4, m12=3
	double mimin(unsigned int i) const;

	//these functions are buggy somewhere!
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

	//!calculate energy of particle \partId in rest frame of system \sys at the invariant mass \invMass_sys
	double eiCms(unsigned int partId, unsigned int sys, double invMass_sys) const;
	//!calculate min value of inv. mass of system \sys given the invariant mass \invMass_sys in system \sys
	double invMassMin(unsigned int sys, unsigned int sys2, double invMass_sys) const;
	//!calculate max value of inv. mass of system \sys given the invariant mass \invMass_sys in system \sys
	double invMassMax(unsigned int sys, unsigned int sys2, double invMass_sys) const;

	//! get mass of paticles
	double getMass(unsigned int num);
	//! get mass of paticles
	double getMass(std::string name);
	//! get spin of decaying particle
	unsigned int getSpin(unsigned int num);
	//! get spin of particles
	unsigned int getSpin(std::string name);

	double getMotherMass() {return M;};
	//! mass of decaying particle
	double M; unsigned int spinM;
	//! width of decaying particle
	double Br;
	//! masses of final state particles
	double m1; unsigned int spin1;
	double m2; unsigned int spin2;
	double m3; unsigned int spin3;
	double m4; unsigned int spin4;
	double m5; unsigned int spin5;
	double m6; unsigned int spin6;

	//! names of particles
	std::string nameMother, name1, name2, name3, name4, name5, name6;

	double m23_sq_min, m23_sq_max;
	double m13_sq_min, m13_sq_max;
	double m12_sq_min, m12_sq_max;
	double m23_min, m23_max;
	double m13_min, m13_max;
	double m12_min, m12_max;

};

#endif
