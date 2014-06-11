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
#include "Core/Kinematics.hpp"
#include "Core/DataPoint.hpp"

void dataPoint::init(){
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	var.reserve(size);
	weight=1;
	return;
}
dataPoint::dataPoint(std::vector<double> vec){
	init();
	var=vec;
	return;
}
dataPoint::dataPoint(Event& ev){
	init();
	Kinematics::instance()->eventToDataPoint(ev,*this);
	return;
}
dataPoint::dataPoint(){
	init();
	return;
}

unsigned int dataPoint::getID(std::string name) const{
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::getVal(): variable with name "<<name<<" not found!";
		return 999;
	}
	return pos;
}

double dataPoint::getVal(std::string name) const{
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::getVal(): variable with name "<<name<<" not found!";
		return -999;
	}
	return getVal(pos);
}
void dataPoint::setVal(std::string name, double val){
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::getVal(): variable with name "<<name<<" not found!";
		return;
	}
	setVal(pos,val);
	return;
}
void dataPoint::setVal(unsigned int num, double val){

	if(var.size()>num+1){
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::setVal(): index for variable out of range!";
		return;
	}
	var[num]=val;
	return;
}
double dataPoint::getVal(unsigned int num) const{

	if(var.size()>num+1){
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::getVal(): index for variable out of range!";
		return 0;
	}
	return var[num];
}
void dataPoint::setPoint(std::vector<double> values){
	unsigned int size = Kinematics::instance()->getVarNames().size();
	if(size!=values.size()){
		BOOST_LOG_TRIVIAL(error)<<"dataPoint2::setPoint(): vector with phsp point out of range!";
		return;
	}
	var=std::vector<double>(values);
	return;
}
std::ostream & operator<<(std::ostream &os, dataPoint &p){
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	os << varNames[0] << "="<<p.getVal(0)<< " | "<<varNames[1]<<"="<<p.getVal(1);
	return os;
}



bool allMasses::Fill(Event &evt){
	dataPoint point(evt);
	// Check number of particle in TClonesrray and if event is within PHSP boundary
	if( nInvMasses != evt.getNParticles()  || !Kinematics::instance()->isWithinPhsp(point))
		return 0;
	eff.push_back(evt.getEfficiency());
	weight.push_back(evt.getWeight());
	sumWeight+=evt.getWeight();

	for(unsigned int pa=0; pa<nInvMasses; pa++){
		for(unsigned int pb=pa+1; pb<nInvMasses; pb++){
			const Particle &inA(evt.getParticle(pa));
			const Particle &inB(evt.getParticle(pb));
			double mymass_sq = inA.invariantMass(inB);
			(masses_sq.at(std::make_pair(pa+1,pb+1))).push_back(mymass_sq);
//			std::cout<<"adding "<<pa+1<<"/"<<pb+1<<" m="<<mymass_sq<<std::endl;
		}//particle loop B
	}//particle loop A
	nEvents++;
	reWeight = (double)nEvents/sumWeight;
	return 1;
}
