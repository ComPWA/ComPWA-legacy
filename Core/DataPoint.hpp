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

class dataPoint
{
private:

public:
	//! Default constructor
	dataPoint();

	/**! Initialize dataPoint with invariant masses.
	 * Missing values are filled by Kinematics
	 */
	dataPoint(int a, int b, double invMassSqA, double invMassSqB);

	//! Construct dataPoint from Event
	dataPoint( const Event& ev );
	//! Construct dataPoint from vector of invariant masses
	dataPoint(std::vector<double> vec);
	~dataPoint(){};

	/**! Fill dataPoint with invariant masses.
	 * Missing values are filled by Kinematics
	 */
	void Set(int a, int b, double invMassSqA, double invMassSqB);
	//! Set value of coordinate name
	void setVal(std::string name, double val);
	//! Set value of coordinate num
	void setVal(unsigned int num, double val);
	//! Get value of coordinate name
	double getVal(std::string name) const;
	//! Get value of coordinate num
	double getVal(unsigned int num) const;

	//! Get ID of coordinate name
	unsigned int getID(std::string name) const;

	//! Set coordinates by vector
	void setPoint(std::vector<double> values);
	//! Get vector of coordinates
	std::vector<double> getPoint() { return var; };

	//! Set weight
	void setWeight(double w) { weight=w; };
	//! Get weight
	double getWeight() { return weight; };
	//! Set efficiency
	void setEfficiency(double e) { eff=e; };
	//! Get efficiency
	double getEfficiency() { return eff; };

	static std::vector<double> getRow(int n, std::vector<dataPoint> v){
		std::vector<double> ret;
		if(!v.size()) return ret;
		if( n >= Kinematics::instance()->GetNVars() )
			throw std::runtime_error("dataPoint::getRow() | out of range!");
		for(int i=0; i<v.size(); i++)
			ret.push_back(v.at(i).getVal(n));
		return ret;
	}

protected:
	void init();
	std::vector<double> var;
	double weight;
	double eff;
	friend std::ostream & operator<<(std::ostream &os, dataPoint &p);
};

#endif /*DPPOINT2_HPP_*/
