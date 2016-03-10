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

#include <algorithm>
#include "Core/Exceptions.hpp"
#include "Core/Kinematics.hpp"
#include "Core/DataPoint.hpp"

void dataPoint::init()
{
	var = std::vector<double>(Kinematics::instance()->GetNVars(), 0);
}

dataPoint::dataPoint(int a, int b, double invMassSqA, double invMassSqB) :
		eff(1.0), weight(1.0)
{
	init();
	Set(a,b,invMassSqA,invMassSqB);
}

void dataPoint::Set(int a, int b, double invMassSqA, double invMassSqB)
{
	Kinematics::instance()->FillDataPoint(a,b,invMassSqA,invMassSqB,*this);
	return;
}

dataPoint::dataPoint(std::vector<double> vec) : weight(1.), eff(1.)
{
	init();
	if(Kinematics::instance()->GetNVars() != vec.size())
		throw std::runtime_error("dataPoint::dataPoint() vector has wrong length!");
	var=vec;
	return;
}

dataPoint::dataPoint( const Event& ev ) : weight(1.), eff(1.)
{
	init();
	Kinematics::instance()->eventToDataPoint(ev,*this);
	weight = ev.getWeight();
	return;
}

dataPoint::dataPoint(): weight(1.), eff(1.)
{
	init();
	return;
}

unsigned int dataPoint::getID(std::string name) const
{
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::getVal() | "
				"Variable with name "<<name<<" not found!";
		return 999;
	}
	return pos;
}

double dataPoint::getVal(std::string name) const
{
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::getVal() | "
				"Variable with name "<<name<<" not found!";
		return -999;
	}
	return getVal(pos);
}

void dataPoint::setVal(std::string name, double val)
{
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::setVal() | "
				"Variable with name "<<name<<" not found!";
		return;
	}
	setVal(pos,val);
	return;
}

void dataPoint::setVal(unsigned int num, double val)
{
	try{
		var.at(num)=val;
	} catch (...) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::setVal() | "
				"Can not access index "<<num<<"!";
		throw;
	}
	return;
}

double dataPoint::getVal(unsigned int num) const
{
	double rt;
	try{
		rt = var.at(num);
	} catch (...) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::getVal() | "
				"Can not access index "<<num<<"!";
		throw;
	}
	return rt;
}

void dataPoint::setPoint(std::vector<double> values)
{
	if(Kinematics::instance()->GetNVars() != values.size())
		throw std::runtime_error("dataPoint::setPoint() vector has wrong length!");
	var=std::vector<double>(values);
	return;
}

std::ostream & operator<<(std::ostream &os, dataPoint &p)
{
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	for(int i=0; i<varNames.size(); i++)
		os << varNames.at(i) << "="<<p.getVal(i)<<" ";
	return os;
}
