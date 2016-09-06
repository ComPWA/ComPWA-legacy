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
#include <stdexcept>

#include "Core/Kinematics.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

void dataPoint::init(){
	var = std::vector<double>(Kinematics::instance()->getVarNames().size(), 0);
}

dataPoint::dataPoint(std::vector<double> vec) : weight(1.){
	init();
	if(Kinematics::instance()->getVarNames().size() != vec.size())
		throw std::runtime_error("dataPoint::dataPoint() vector has wrong length!");
	var=vec;
	return;
}
dataPoint::dataPoint(const Event& ev): weight(1.){
	init();
	Kinematics::instance()->eventToDataPoint(ev,*this);
	return;
}
dataPoint::dataPoint(): weight(1.){
	init();
	return;
}

void dataPoint::reset(unsigned int size) {
  var = std::vector<double>(size);
}

unsigned int dataPoint::getID(std::string name) const{
  return Kinematics::instance()->getVariableIndex(name);
}

double dataPoint::getVal(std::string name) const{
	return getVal(getID(name));
}
void dataPoint::setVal(std::string name, double val){
	setVal(getID(name),val);
}
void dataPoint::setVal(unsigned int num, double val){
	try{
		var.at(num)=val;
	} catch (...) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::setVal(): cannot access index "<<num<<"!";
		throw;
	}
	return;
}
double dataPoint::getVal(unsigned int num) const{
  return var[num];
	/*double rt;
	try{
		rt = var.at(num);
	} catch (...) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::getVal(): cannot access index "<<num<<"!";
		throw;
	}
	return rt;*/
}
void dataPoint::setPoint(std::vector<double> values){
	if(Kinematics::instance()->getVarNames().size() != values.size())
		throw std::runtime_error("dataPoint::setPoint() vector has wrong length!");
	var=std::vector<double>(values);
	return;
}
std::ostream & operator<<(std::ostream &os, dataPoint &p){
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	for(int i=0; i<varNames.size(); i++)
		os << varNames.at(i) << "="<<p.getVal(i)<<" ";
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

} /* namespace ComPWA */
