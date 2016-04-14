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
#include <memory>

#include "Core/Logging.hpp"

#include "Core/Event.hpp"
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
	virtual void EventToDataPoint(const Event& ev, dataPoint& point) const = 0;

	//! Event to dataPoint conversion
	virtual void FillDataPoint(int a, int b, double invMassSqA, double invMassSqB,
			dataPoint& point) = 0;

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
	virtual bool IsWithinPhsp(const dataPoint& point) { };

	/**! Checks if the position is within the phase-space boundaries.
	 * This only works correctly if both variables are orthogonal to each other.
	 * E.g. and invariant mass and an angle.
	 * @param idA Variable id of invariant mass
	 * @param idB Variable id of angle
	 * @param varA Invariant mass
	 * @param varB Helicity angle
	 * @return
	 */
	virtual bool IsWithinBoxPhsp(int idA, int idB, double varA, double varB) = 0;

	//! Get name of mother particle
	virtual std::string GetMotherName() { return _nameMother; };
	//! Get mass of mother particle
	virtual double GetMotherMass() { return _M; }
	//! calculated the PHSP volume of the current decay by MC integration
	virtual double GetPhspVolume() = 0;


	//! get mass of particles
	virtual double GetMass(unsigned int num) = 0;
	//! get mass of paticles
	virtual double GetMass(std::string name) = 0;
	//! Get number of particles
	virtual unsigned int GetNumberOfParticles() { return _nPart; }
	//! Get number of variables
	virtual unsigned int GetNVars() const { return _varNames.size(); }

	/** Calculate Break-up momentum squared
	 *
	 * Calculate Break-up momentum at energy @sqrtS for particles with masses @ma and @mb.
	 * From PDG2014 Eq.46-20a. Below threshold the function is analytically continued.
	 * @param sqrtS center-of-mass energy
	 * @param ma mass particle A
	 * @param mb mass particle B
	 * @return |break-up momentum|
	 */
	static double qSqValue(double sqrtS, double ma, double mb);
	/** Calculate Break-up momentum
	 *
	 * Calculate Break-up momentum at energy @sqrtS for particles with masses @ma and @mb.
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
		_nameMother(nameM), _Br(widthM), _nPart(n) {};

	//!Copy constructor (protected)
	Kinematics(const Kinematics&) {};

	//! Default destructor
	virtual ~Kinematics() {};

};

class TwoBodyKinematics : public Kinematics
{
public:
	TwoBodyKinematics(std::string _nameMother, std::string _name1,
			std::string _name2, double deltaMassWindow=0.5);

	void init();

	static Kinematics* createInstance(std::string _nameMother,
			std::string _name1, std::string _name2, double massWindow=0.5){
		if(_inst) return _inst;
		_inst = new TwoBodyKinematics(_nameMother, _name1, _name2, massWindow);
		return _inst;
	}

	//! Converts Event to dataPoint
	virtual void EventToDataPoint(const Event& ev, dataPoint& point) const;

	virtual void FillDataPoint(int a, int b, double invMassSqA, double invMassSqB,
			dataPoint& point) { };

	//! checks of data point is within phase space boundaries
	virtual bool IsWithinPhsp(const dataPoint& point);

	/**! Checks if the position is within the phase-space boundaries.
	 * This only works correctly if both variables are orthogonal to each other.
	 * E.g. and invariant mass and an angle.
	 * @param idA Variable id of invariant mass
	 * @param idB Variable id of angle
	 * @param varA Invariant mass
	 * @param varB Helicity angle
	 * @return
	 */
	virtual bool IsWithinBoxPhsp(int idA, int idB, double varA, double varB) { };

	//! Calculate phase-space volume
	virtual double GetPhspVolume() { return (mass_max-mass_min); }

	//! get mass of particles
	virtual double GetMass(unsigned int num);

	//! get mass of paticles
	virtual double GetMass(std::string name);

protected:
	std::string name1;//! name of daughter 1
	double mSq1; //! masse squared of daughter 1
	double m1; //! masses of daughter 1
	unsigned int spin1; //! spin of daughter 1

	std::string name2;//! name of daughter 2
	double mSq2; //! masse squared of daughter 2
	double m2; //! masses of daughter 2
	unsigned int spin2;//! spin of daughter 2

	double mass_sq_min; //!minimum value of masssq
	double mass_sq_max;//!maximum value of masssq
	double mass_min; //!minimum value of masssq
	double mass_max;//!maximum value of masssq

};
#endif /* KINEMATICS_HPP_ */
