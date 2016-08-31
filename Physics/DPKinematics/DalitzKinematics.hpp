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

#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace DPKinematics {

class DalitzKinematics : public Kinematics
{
protected:
	bool massIdsSet;
	unsigned int id23;
	unsigned int id13;

	//! calculated dalitz plot area for the given kinematics
	double calculatePSArea();
	//! initialization
	void init();
	//! default constructor
	DalitzKinematics():massIdsSet(false){};
	//! constructor access particles by name, masses etc are obtained from PhysConst singleton
	DalitzKinematics(std::string _nameMother, std::string _name1, std::string _name2, std::string _name3);
	//! constructor with particle masses and names, independent from PhysConst
	DalitzKinematics(double _M, double _Br, double _m1, double _m2, double _m3,
			std::string _nameMother, std::string _name1, std::string _name2, std::string _name3);

  //! Event to dataPoint conversion
  void translateEventToDataPoint(const Event& ev, dataPoint& point) const;

public:
	static Kinematics* createInstance(std::string _nameMother, std::string _name1, std::string _name2, std::string _name3){
		if(0 == inst_)
	    inst_ = new DalitzKinematics(_nameMother, _name1, _name2,_name3);
		return inst_;
	}

	/**
	 * \brief Generate contour of phsp boundary
	 *
	 * @param xsys Which subsystem should be plotted on x?
	 * @param ysys Which subsystem should be plotted on y?
	 * @param n Number of points
	 * @param xcoord array with x values
	 * @param ycoord array with y values
	 *
	 * The allocated size of the arrays should be n+1.
	 */
	void phspContour(unsigned int xsys,unsigned int ysys, unsigned int n, double* xcoord, double* ycoord);
	/*! Calculates the helicity angle.
	 *
	 * Calculates the helicity angle for subsystem @param sys given the invariant masses
	 * @param invMass23sq and @param invMass23sq. The angle is measured versus daughter 2 in system [12],
	 * versus daughter 1 in [13] and versus 2 in [23]
	 */
	double helicityAngle(unsigned int sys, double invMassSq23, double invMassSq13);
	//! Helicity angle for subSystem sys at dataPoint point.
	double helicityAngle(unsigned int sys, const dataPoint& point);
	/**
	 * \brief Calculates the scattering angle.
	 *
	 * Function obsolete!
	 * Calculates the scattering angle given the invariant masses @param s and @param t.
	 * The angle is measured between the spectator particle @param mSpec and particle @param m.
	 * @param mSecond is the third particle of the decay
	 * You should set the masses as follows:
	 * (m12sq,m23sq,M,m3,m1,m2); for subsystem 3
	 * (m13sq,m12sq,M,m2,m3,m1); for subsystem 4
	 * (m23sq,m13sq,M,m1,m2,m3); for subsystem 5
	 * When masses for scatteringAngle() are set correctly both functions are equivalent.
	 * \return cos(helicityAngle)
	 */
	double scatteringAngle(double s, double t, double M, double mSpec, double mSecond, double m);
	//! Calculates third dalitz plot variable, e.g f(s1,s2)=s3
	double getThirdVariableSq(double, double) const;
	//! Checks if data point is within phase space boundaries
	bool isWithinPhsp(const dataPoint &point) ;
	//! Returns the dalitz plot area for the given kinematics
	//double getPhspVolume();
    //! Returns the dalitz plot area for the given kinematics and limited m23 range
    double getPhspVolumePart(double, double);
	//! Calculated momenta n,m using legendre polynomials
	double calculateMoments(unsigned int sys, dataPoint& point, unsigned int n, unsigned int m);
	//!maximum value for invariant mass squared: m23sq=5, m13sq=4, m12sq=3
	double mimax(unsigned int sys) const;
	//!minimum value for invariant mass squared: m23sq=5, m13sq=4, m12sq=3
	double mimin(unsigned int sys) const;

	//these functions are buggy somewhere!
	//	double lambda(double x, double y, double z)const;
	//	double s2min(double s1, double m0, double m1, double m2, double m3)const;
	//	double s2max(double s1, double m0, double m1, double m2, double m3)const;
	//	double s3min(double s1, double m0, double m1, double m2, double m3)const;
	//	double s3max(double s1, double m0, double m1, double m2, double m3)const;
	//	double s1min(double s2, double m0, double m1, double m2, double m3)const;
	//	double s1max(double s2, double m0, double m1, double m2, double m3)const;
	//	double s2min(double s1)const { return s2min(s1,M,m1,m2,m3); };
	//	double s2max(double s1)const { return s2max(s1,M,m1,m2,m3); };
	//	double s3min(double s1)const { return s3min(s1,M,m1,m2,m3); };
	//	double s3max(double s1)const { return s3max(s1,M,m1,m2,m3); };
	//	double s1min(double s2)const { return s1min(s2,M,m1,m2,m3); };
	//	double s1max(double s2)const { return s1max(s2,M,m1,m2,m3); };

	//!calculate energy of particle partId in rest frame of system sys at the invariant mass invMass_sys
	double eiCms(unsigned int partId, unsigned int sys, double invMass_sys) const;
	//!calculate min value of inv. mass of system sys2 given the invariant mass invMass_sys in system sys
	double invMassMin(unsigned int sys, unsigned int sys2, double invMass_sys) const;
	//!calculate max value of inv. mass of system sys2 given the invariant mass invMass_sys in system sys
	double invMassMax(unsigned int sys, unsigned int sys2, double invMass_sys) const;

	//! returns absolute minimum for variable
	double getMin(std::string);
	//! returns absolute maximum for variable
	double getMax(std::string);
	//! get mass of paticles
	double getMass(unsigned int num) const;
	//! get mass of paticles
	double getMass(std::string name) const;
	//! get spin of decaying particle
	unsigned int getSpin(unsigned int num);
	//! get spin of particles
	unsigned int getSpin(std::string name);

	//! mass of mother particle
	double getMotherMass() const {return M;};

	std::string nameMother;//! name of mother particle
	double Msq; //! mass squared of mother particle
	double M; //! mass of mother particle
	unsigned int spinM;//! spin of mother particle
	double Br;//! width of decaying particle

	std::string name1;//! name of daughter 1
	double mSq1; //! masse squared of daughter 1
	double m1; //! masses of daughter 1
	unsigned int spin1; //! spin of daughter 1
	std::string name2;//! name of daughter 2
	double mSq2; //! masse squared of daughter 2
	double m2; //! masses of daughter 2
	unsigned int spin2;//! spin of daughter 2
	std::string name3;//! name of daughter 3
	double mSq3; //! masse squared of daughter 3
	double m3; //! masses of daughter 3
	unsigned int spin3;//! spin of daughter 3
	std::string name4;
	double mSq4; //! masse squared of daughter 4
	double m4;
	unsigned int spin4;


	double m23_sq_min; //!minimum value of m23sq
	double m23_sq_max;//!maximum value of m23sq
	double m13_sq_min;//!minimum value of m13sq
	double m13_sq_max;//!maximum value of m13sq
	double m12_sq_min;//!minimum value of m12sq
	double m12_sq_max;//!maximum value of m12sq
	double m23_min; //!minimum value of m23sq
	double m23_max;//!maximum value of m23sq
	double m13_min; //!minimum value of m13sq
	double m13_max;//!maximum value of m13sq
	double m12_min; //!minimum value of m12sq
	double m12_max;//!maximum value of m12sq

};

} /* namespace DPKinematics */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
