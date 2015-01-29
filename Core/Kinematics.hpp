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


protected:
	std::vector<std::string> varNames;
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
	static Kinematics* createInstance(std::string _nameMother, std::string _name1, std::string _name2){
		_inst = new TwoBodyKinematics(_nameMother, _name1, _name2);
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
