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
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Physics/DPKinematics/DataPoint2.hpp"

dataPoint2::dataPoint2(){
	std::vector<std::string> varNames = DalitzKinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	var.reserve(size);
	return;

}

double dataPoint2::getVal(std::string name){
	std::vector<std::string> varNames = DalitzKinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::getVal(): variable with name "<<name<<" not found!";
		return -999;
	}
	return getVal(pos);
}
void dataPoint2::setVal(std::string name, double val){
	std::vector<std::string> varNames = DalitzKinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::getVal(): variable with name "<<name<<" not found!";
		return;
	}
	setVal(pos,val);
	return;
}
void dataPoint2::setVal(unsigned int num, double val){

	if(var.size()>num+1){
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::setVal(): index for variable out of range!";
		return;
	}
	var[num]=val;
	return;
}
double dataPoint2::getVal(unsigned int num){

	if(var.size()>num+1){
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::getVal(): index for variable out of range!";
		return 0;
	}
	return var[num];
}
void dataPoint2::setPoint(std::vector<double> values){
	unsigned int size = DalitzKinematics::instance()->getVarNames().size();
	if(size!=values.size()){
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::setPoint(): vector with phsp point out of range!";
		return;
	}
	var=std::vector<double>(values);
	return;
}
