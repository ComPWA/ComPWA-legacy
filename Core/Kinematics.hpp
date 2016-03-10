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

static const char * formFactorTypeString[] = { "noFormFactor", "BlattWeisskopf", "CrystalBarrel" };
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
	//! vector with names of variables, e.g. vec[0]=m23sq, vec[1]=m13sq
	const std::vector<std::string>& getVarNames(){return varNames;}
	//! checks of data point is within phase space boundaries
	virtual bool isWithinPhsp(const dataPoint& point) = 0;
	//! mass of mother particle
	virtual double getMotherMass() = 0;
	//! calculated the PHSP volume of the current decay by MC integration
	virtual double getPhspVolume() = 0;
	//! converts Event to dataPoint
	virtual void eventToDataPoint(Event& ev, dataPoint& point) = 0;
	//! Event to dataPoint conversion
	virtual void FillDataPoint(int a, int b, double invMassSqA, double invMassSqB,
			dataPoint& point) = 0;

	//! get mass of particles
	virtual double getMass(unsigned int num) = 0;
	//! get mass of paticles
	virtual double getMass(std::string name) = 0;
	//! Get number of particles
	virtual unsigned int getNumberOfParticles() { return nPart; }
	//! Get number of variables
	virtual unsigned int GetNVars() const { return varNames.size(); }

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
	static double FormFactor(double sqrtS, double ma, double mb, double spin, double mesonRadius,
			formFactorType type=formFactorType::BlattWeisskopf);

protected:
	unsigned int nPart;
	//! Internal names of variabes
	std::vector<std::string> varNames;
	//! Latex titles for variables
	std::vector<std::string> varTitles;
	static Kinematics* _inst;
	Kinematics() {};
	Kinematics(const Kinematics&) {};
	virtual ~Kinematics() {};
};

class TwoBodyKinematics : public Kinematics
{
public:
	TwoBodyKinematics(std::string _nameMother, std::string _name1, std::string _name2, double deltaMassWindow=0.5);
	void init();
	static Kinematics* createInstance(std::string _nameMother, std::string _name1, std::string _name2, double massWindow=0.5){
		if(_inst) return _inst;
		_inst = new TwoBodyKinematics(_nameMother, _name1, _name2, massWindow);
		return _inst;
	}
	//! checks of data point is within phase space boundaries
	virtual bool isWithinPhsp(const dataPoint& point);
	//! mass of mother particle
	virtual double getMotherMass() { return M; }
	//! calculated the PHSP volume of the current decay by MC integration
	virtual double getPhspVolume() { return (mass_max-mass_min); }
	//! converts Event to dataPoint
	virtual void eventToDataPoint(Event& ev, dataPoint& point);
	virtual void FillDataPoint(int a, int b, double invMassSqA, double invMassSqB,
			dataPoint& point) { };
	//! get mass of particles
	virtual double getMass(unsigned int num);
	//! get mass of paticles
	virtual double getMass(std::string name);

protected:
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

	double mass_sq_min; //!minimum value of masssq
	double mass_sq_max;//!maximum value of masssq
	double mass_min; //!minimum value of masssq
	double mass_max;//!maximum value of masssq

};
#endif /* KINEMATICS_HPP_ */
