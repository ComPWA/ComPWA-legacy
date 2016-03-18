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

public:
	static Kinematics* createInstance(std::string nameMother,
			std::string name1, std::string name2, std::string name3)
	{
		if(_inst) return _inst;
		_inst = new DalitzKinematics(nameMother, name1, name2, name3);
		return _inst;
	}

	//! Event to dataPoint conversion
	void EventToDataPoint(const Event& ev, dataPoint& point) const;

	//! Event to dataPoint conversion
	void FillDataPoint(int a, int b, double invMassSqA, double invMassSqB,
			dataPoint& point);

	/**! Generate contour of phsp boundary
	 *
	 * @param xsys Which subsystem should be plotted on x?
	 * @param ysys Which subsystem should be plotted on y?
	 * @param n Number of points
	 * @param xcoord array with x values
	 * @param ycoord array with y values
	 *
	 * The allocated size of the arrays should be n+1.
	 */
	void phspContour(unsigned int xsys,unsigned int ysys, unsigned int n,
			double* xcoord, double* ycoord);

	/**! Calculates the helicity angle.
	 *
	 * Calculates the helicity angle for subsystem @param sys given
	 * the invariant masses. @param invMass23sq and @param invMass23sq. The
	 * angle is measured versus daughter 2 in system [12],
	 * versus daughter 1 in [13] and versus 2 in [23]
	 */
	double helicityAngle(unsigned int sys, double invMassSq23, double invMassSq13);

	/**! Calculates the helicity angle of particle m and mSpec in the rest
	 * Helicity angle between particle m and mSpec in the rest frame of m and m2
	 * @param M Mass of mother particle
	 * @param m Mass of particle
	 * @param m2 Mass of recoiling particle
	 * @param mSpec Mass of spectator particle
	 * @param invMassSqA invariant mass of m+m2
	 * @param invMassSqB invariant mass of m+mSpec
	 * @return
	 */
	static double helicityAngle(double M, double m, double m2, double mSpec,
		double invMassSqA, double invMassSqB);

	//! Helicity angle for subSystem sys at dataPoint point.
	double helicityAngle(unsigned int sys,dataPoint& point);

	/**! Calculates the scattering angle.
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
	bool IsWithinPhsp(const dataPoint &point) ;

	/**! Checks if the position is within the phase-space boundaries.
	 * This only works correctly if both variables are orthogonal to each other.
	 * E.g. and invariant mass and an angle.
	 * @param idA Variable id of invariant mass
	 * @param idB Variable id of angle
	 * @param varA Invariant mass
	 * @param varB Helicity angle
	 * @return
	 */
	bool IsWithinBoxPhsp(int idA, int idB, double varA, double varB);

	//! Returns the dalitz plot area for the given kinematics
	double GetPhspVolume();

	//! Calculated momenta n,m using legendre polynomials
	double calculateMoments(unsigned int sys, dataPoint& point, unsigned int n, unsigned int m);

	//! Global minimum and maximum value of variable
	std::pair<double,double> GetMinMax(std::string varName) const;

	//! Global minimum and maximum value of variable
	std::pair<double,double> GetMinMax(unsigned int sys) const;

	//! Energy of particle partId in rest frame of system sys at the invariant mass invMass_sys
	double eiCms(unsigned int partId, unsigned int sys, double invMass_sys) const;

	//! Local maximum and minimum value of variable
	std::pair<double,double> GetMinMaxLocal(unsigned int sys, unsigned int sys2,
			double invMass_sys) const;

	//! get mass of paticles
	double GetMass(unsigned int num);
	//! get mass of paticles
	double GetMass(std::string name);
	//! get spin of decaying particle
	unsigned int getSpin(unsigned int num);
	//! get spin of particles
	unsigned int getSpin(std::string name);

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

protected:
	//! default constructor
	DalitzKinematics() : massIdsSet(false) { };

	//! constructor access particles by name, masses etc are obtained from PhysConst singleton
	DalitzKinematics(std::string _nameMother, std::string _name1,
			std::string _name2, std::string _name3);

	//! constructor with particle masses and names, independent from PhysConst
	DalitzKinematics(double _M, double _Br, double _m1, double _m2, double _m3,
			std::string _nameMother, std::string _name1,
			std::string _name2, std::string _name3);

	//! initialization
	void init();

	// Check if variables are orthogonal to each other
	bool AreBoxVariables(unsigned int idA, unsigned int idB);

	//! calculated dalitz plot area for the given kinematics
	void calcDParea();

	bool massIdsSet;
	bool _DPareaCalculated;	//! is phsp area already calculated?
	double _DParea;	//! phsp area

};

#endif
