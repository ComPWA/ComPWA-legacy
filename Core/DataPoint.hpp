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
#include "Core/Kinematics.hpp"
class Kinematics;

class dataPoint
{
private:

public:

	dataPoint(dataPoint const&){};
	dataPoint();
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
	void setVal(unsigned int num, double val);
	double getVal(unsigned int num) const;
	void setPoint(std::vector<double> values);

protected:
	std::vector<double> var;
};


#endif /* DPPOINT2_HPP_ */
