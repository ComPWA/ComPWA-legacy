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
#ifndef DALITZEFFICIENCY_HPP_
#define DALITZEFFICIENCY_HPP_

#include <iostream>
#include <vector>
#include <boost/log/trivial.hpp>
using namespace boost::log;

/**
 *  \class DalitzEfficiency
 *  \brief Virtual efficiency class
 */
class Efficiency{
private:

public:
	Efficiency();

//	DalitzEfficiency(const DalitzEfficiency&, const char*);

	virtual ~Efficiency();

	virtual double evaluate(std::vector<double> x) = 0;
};

/**
 *  \class UnitEfficiency
 *  \brief implementation of virtual class efficiency. Efficiency ist constant one allover the PHSP
 */
class UnitEfficiency : public Efficiency {
private:
public:
	UnitEfficiency(){
		BOOST_LOG_TRIVIAL(info)<<"Efficiency: creating UnitEfficiency!";
	};
	~UnitEfficiency(){};
	virtual double evaluate(std::vector<double> x) {return 1;};
};
#endif /* DALITZEFFICIENCY_HPP_ */
