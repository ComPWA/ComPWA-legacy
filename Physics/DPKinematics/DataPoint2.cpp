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
	varNames = DalitzKinematics::instance()->getVarNames();
	size=varNames.size();
	var.reserve(size);
	return;

}
dataPoint2::dataPoint2(std::vector<std::string> names): varNames(names){
	size=varNames.size();
	var.reserve(size);
	return;

}
double dataPoint2::getVal(std::string name){
	int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		std::cout<<"ERROR: dataPoint2: getVal(): variable with name"<<name<<" not found!"<<std::endl;
		return -999;
	}
	return getVal(pos);
}
void dataPoint2::setVal(std::string name, double val){
	int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		std::cout<<"ERROR: dataPoint2: setVal(): variable with name"<<name<<" not found!"<<std::endl;
		return;
	}
	setVal(pos,val);
	return;
}
void dataPoint2::setVal(unsigned int num, double val){

	if(var.size()>num+1){
		std::cout<<"ERROR: dataPoint2: setVal(): index for variable out of range!"<<std::endl;
		return;
	}
	var[num]=val;
	return;
}
double dataPoint2::getVal(unsigned int num){

	if(var.size()>num+1){
		std::cout<<"ERROR: dataPoint2: getVal(): index for variable out of range!"<<std::endl;
		return 0;
	}
	return var[num];
}
void dataPoint2::setPoint(std::vector<double> values){
	if(varNames.size()!=values.size()){
		std::cout<<"ERROR: dataPoint2: setPoint(): vector with phsp point out of range!"<<std::endl;
		return;
	}
	var=std::vector<double>(values);
	return;
}
