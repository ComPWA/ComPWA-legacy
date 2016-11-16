//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff -
//-------------------------------------------------------------------------------

#ifndef KINEMATICS_HPP_
#define KINEMATICS_HPP_

#include <vector>
#include <string>
#include <complex>

#include <Core/Event.hpp>

namespace ComPWA {

class dataPoint;

static const char * formFactorTypeString[] = {
		"noFormFactor",
		"BlattWeisskopf",
		"CrystalBarrel"
};

enum formFactorType{
	noFormFactor = 0,
	BlattWeisskopf = 1,
	CrystalBarrel = 2
};

class Kinematics
{
public:
	//! singleton pattern
	static Kinematics* instance();

	//! converts Event to dataPoint
	virtual void EventToDataPoint( const ComPWA::Event& ev,
			dataPoint& point) const = 0;

	//! Event to dataPoint conversion
	virtual void FillDataPoint(int a, int b, double invMassSqA,
			double invMassSqB, dataPoint& point) const = 0;

	//! vector with names of variables, e.g. vec[0]=m23sq, vec[1]=m13sq
	std::vector<std::string> GetVarNames() const { return _varNames; }

	//! vector with names of variables, e.g. vec[0]=m23sq, vec[1]=m13sq
	std::string GetVarName(unsigned int pos) const { return _varNames.at(pos); }

	//! vector with names of variables, e.g. vec[0]=m23sq, vec[1]=m13sq
	std::vector<std::string> GetVarTitles() const { return _varTitles; }

	//! vector with names of variables, e.g. vec[0]=m23sq, vec[1]=m13sq
	std::string GetVarTitle(unsigned int pos) const { return _varTitles.at(pos); }

	//! Get position of variable @param varName
	unsigned int FindVariable(std::string varName) const;

	//! Checks of data point is within phase space boundaries
	virtual bool IsWithinPhsp(const dataPoint& point) const { };

	/**! Checks if the position is within the phase-space boundaries.
	 * This only works correctly if both variables are orthogonal to each other.
	 * E.g. and invariant mass and an angle.
	 * @param idA Variable id of invariant mass
	 * @param idB Variable id of angle
	 * @param varA Invariant mass
	 * @param varB Helicity angle
	 * @return
	 */
	virtual bool IsWithinBoxPhsp(int idA, int idB, double varA, double varB) const = 0;

	//! Get name of mother particle
	virtual std::string GetMotherName() const { return _nameMother; };
	//! Get mass of mother particle
	virtual double GetMotherMass() const { return _M; }
	//! calculated the PHSP volume of the current decay by MC integration
	virtual double GetPhspVolume();

	//! get mass of particles
	virtual double GetMass(unsigned int num) const = 0;
	//! get mass of particles
	virtual double GetMass(std::string name) const = 0;
	//! Get number of particles
	virtual unsigned int GetNumberOfParticles() const { return _nPart; }
	//! Get number of variables
	virtual unsigned int GetNVars() const { return _varNames.size(); }

	/** Calculate Break-up momentum squared
	 *
	 * Calculate Break-up momentum at energy @param sqrtS for particles with masses @param ma and @param mb .
	 * From PDG2014 Eq.46-20a. Below threshold the function is analytically continued.
	 * @param sqrtS center-of-mass energy
	 * @param ma mass particle A
	 * @param mb mass particle B
	 * @return |break-up momentum|
	 */
	static double qSqValue(double sqrtS, double ma, double mb);
	/** Calculate Break-up momentum
	 *
	 * Calculate Break-up momentum at energy @param sqrtS for particles with masses @param ma and @param mb .
	 * From PDG2014 Eq.46-20a. Below threshold the function is analytically continued.
	 * @param sqrtS center-of-mass energy
	 * @param ma mass particle A
	 * @param mb mass particle B
	 * @return |break-up momentum|
	 */
	static std::complex<double> qValue(double sqrtS, double ma, double mb);
	/** Two body phsp factor
	 *
	 * From PDG2014 Eqn.47-2
	 * @param sqrtS invariant mass of particles A and B
	 * @param ma Mass of particle A
	 * @param mb Mass of particle B
	 * @return
	 */
	static std::complex<double> phspFactor(double sqrtS, double ma, double mb);

	//! Calculate form factor
	static double FormFactor(double sqrtS, double ma, double mb,
			double spin, double mesonRadius,
			formFactorType type=formFactorType::BlattWeisskopf);

	//! Calculate form factor
	static double FormFactor(double sqrtS, double ma, double mb,
			double spin, double mesonRadius, std::complex<double> qValue,
			formFactorType type=formFactorType::BlattWeisskopf);

protected:
	//! Number of particles in reaction
	unsigned int _nPart;

	//! Internal names of variabes
	std::vector<std::string> _varNames;
	//! Latex titles for variables
	std::vector<std::string> _varTitles;

	//Parameters of decaying mother particle (we assume that we have a decay)
	std::string _nameMother;//! name of mother particle
	double _M; //! mass of mother particle
	double _Msq; //! mass of mother particle
	unsigned int _spinM;//! spin of mother particle
	double _Br;//! width of decaying particle

	//Singleton stuff
	static Kinematics* _inst;

	//!Default constructor (protected)
	Kinematics(std::string nameM="", double widthM=0.0, unsigned int n=3) :
		_nameMother(nameM), _Br(widthM), _nPart(n),
		is_PS_area_calculated_(false), PS_area_(0.0),
		_M(1.0), _Msq(1.0), _spinM(0) { };

	//!Delete Copy constructor (protected)
	Kinematics(const Kinematics&) = delete;

	//! Default destructor
	virtual ~Kinematics() {};

	//! Delete assignment operator
	void operator=(const Kinematics&) = delete;

	virtual double calculatePSArea() =0;
	bool is_PS_area_calculated_;
	double PS_area_;
};

} /* namespace ComPWA */
#endif /* KINEMATICS_HPP_ */
