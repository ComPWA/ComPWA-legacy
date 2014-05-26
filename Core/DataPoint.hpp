//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

//! dataPoint stores kinematic information of a dalitz plot
/*!
 * @file DataPoint.hpp
 *\class dataPoint
 *      dataPoint is a singleton class which provides a
 *      certain phase space point to all classes of the framework. The point can be set anywhere in
 *      the framework and can be read by any amplitude class.
 */

#ifndef DPPOINT2_HPP_
#define DPPOINT2_HPP_

#include <cstdlib>
#include <math.h>
#include <vector>
#include <map>
#include "Core/Kinematics.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Event.hpp"

class allMasses
{
public:
	/**Constructor
	 *
	 * @param inMasses number of variables 3 for 3-body decay
	 * @param inTup combinations for inv masses (e.g (2,3), (1,3), (1,2))
	 */
	allMasses(unsigned int inMasses, std::vector<std::pair<unsigned int, unsigned int> >& inTup) :
		nInvMasses(inMasses),nEvents(0) {
		for(unsigned int i=0; i<inTup.size(); i++)
			masses_sq.insert( std::make_pair( inTup[i], std::vector<double>() ) );
	}


	/**Constructor
	 * additionally allocates memory, if number of events is known.
	 *
	 * @param inMasses number of variables 3 for 3-body decay
	 * @param inEvents number of events to reserve memory for
	 * @param inTup combinations for inv masses (e.g (2,3), (1,3), (1,2))
	 */
	allMasses(unsigned int inMasses, unsigned int inEvents, std::vector<std::pair<unsigned int, unsigned int> >& inTup)
	:nInvMasses(inMasses),nEvents(inEvents){
		//alocate memory in advance
		for(unsigned int i=0; i<inTup.size(); i++)
			masses_sq.insert( std::make_pair( inTup[i], std::vector<double>(inEvents,0.) ) );
		eff = std::vector<double>(nEvents,1.);
		weight = std::vector<double>(nEvents,1.);
	}

	//!Standard constructor
	allMasses():nInvMasses(0),nEvents(0) {}

	//! Fill event
	bool Fill(Event &evt);
	void setEfficiency(std::shared_ptr<Efficiency> effObj){
		unsigned int nEvents = masses_sq.size();
		eff = std::vector<double>(nEvents,1.);
		for(unsigned int i=0; i<nEvents;i++){
			std::vector<double> data;
			data.push_back(masses_sq.at( std::make_pair(2,3) )[i]);
			data.push_back(masses_sq.at( std::make_pair(1,3) )[i]);
			double value  = effObj->evaluate(data);
			//			if(value <= 0) value = 0;
			//			else value = 1/value;
			/*
			 * we need to use sqrt(eff) here because in the current
			 * implementation the Amplitude value is squared after
			 * multiplication with the efficiency
			 */
			eff.at(i) = sqrt(value);
			if(value==0) eff.at(i) = 0.001;
			//			std::cout<<effObj->evaluate(data)<<std::endl;
		}
		return;
	}

	std::map<std::pair<unsigned int, unsigned int>,std::vector<double> > masses_sq;
	std::vector<double> eff;
	std::vector<double> weight;
	unsigned int nEvents;
	unsigned int nInvMasses;

};

class dataPoint
{
private:

public:

	dataPoint(dataPoint const&){};
	dataPoint();
	dataPoint(Event& ev);
	dataPoint(std::vector<double> vec);
	~dataPoint(){};
	//! checks if point lies within phase space boundaries
	//	bool isWithinPhsp() const{ return DalitzKinematics::instance()->isWithinPhsp(this); };
	//
	//	//! get inv. mass for subsys: 1+2=3, 1+3=4, 2+3=5
	//	double getM(int subsys) {return sqrt(getMsq(subsys));};
	//	//! get inv. mass of particle a and b
	//	double getM(int a, int b) {return getM(a+b);};//daughter1 and daughter2
	//	//! get inv. mass sq. for subsys: 1+2=3, 1+3=4, 2+3=5
	//	double getMsq(int subsys);
	//	//! get inv. mass sq. of particle a and b
	//	double getMsq(int a, int b) {return getMsq(a+b);};//daughter1 and daughter2
	//	//! set inv. mass for subsys: 1+2=3, 1+3=4, 2+3=5
	//	void setM(int sys, double val) { setMsq(sys, val*val); };
	//	//! set inv. mass of particle a and b
	//	void setM(int a, int b, double val) { setM(a+b, val);};
	//	//! set inv. mass sq. for subsys: 1+2=3, 1+3=4, 2+3=5
	//	void setMsq(int, double);
	//	//! set inv. mass sq. of particle a and b
	//	void setMsq(int a, int b, double val) { setMsq(a+b, val);};

	void setVal(std::string name, double val);
	double getVal(std::string name) const;
	unsigned int getID(std::string name) const;
	void setVal(unsigned int num, double val);
	double getVal(unsigned int num) const;
	void setPoint(std::vector<double> values);
	void setWeight(double w) { weight=w; };
	double getWeight() { return weight; };

protected:
	void init();
	std::vector<double> var;
	double weight;
	friend std::ostream & operator<<(std::ostream &os, dataPoint &p);
};
#endif /*DPPOINT2_HPP_*/
