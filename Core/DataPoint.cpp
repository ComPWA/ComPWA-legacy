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
	var = std::vector<double>(Kinematics::instance()->getVarNames().size(), 0);
}

dataPoint::dataPoint(std::vector<double> vec) : weight(1.){
	init();
	if(Kinematics::instance()->getVarNames().size() != vec.size())
		throw std::runtime_error("dataPoint::dataPoint() vector has wrong length!");
	var=vec;
	return;
}
dataPoint::dataPoint(Event& ev): weight(1.){
	init();
	Kinematics::instance()->eventToDataPoint(ev,*this);
	return;
}
dataPoint::dataPoint(): weight(1.){
	init();
	return;
}

unsigned int dataPoint::getID(std::string name) const{
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::getVal(): variable with name "<<name<<" not found!";
		return 999;
	}
	return pos;
}

double dataPoint::getVal(std::string name) const{
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::getVal(): variable with name "<<name<<" not found!";
		return -999;
	}
	return getVal(pos);
}
void dataPoint::setVal(std::string name, double val){
	std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
	unsigned int size = varNames.size();
	unsigned int pos = find(varNames.begin(), varNames.end(), name) - varNames.begin();
	if(pos<0||pos>size-1) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::setVal(): variable with name "<<name<<" not found!";
		return;
	}
	setVal(pos,val);
	return;
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
	double rt;
	try{
		rt = var.at(num);
	} catch (...) {
		BOOST_LOG_TRIVIAL(error)<<"dataPoint::getVal(): cannot access index "<<num<<"!";
		throw;
	}
	return rt;
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


//--------------------------------------------------------

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
		}//particle loop B
	}//particle loop A
	nEvents++;
	reWeight = (double)nEvents/sumWeight;
	return 1;
}


allMasses::allMasses(unsigned int inMasses, std::vector<std::pair<unsigned int, unsigned int> >& inTup) :
				nInvMasses(inMasses),nEvents(0) {
	for(unsigned int i=0; i<inTup.size(); i++)
		masses_sq.insert( std::make_pair( inTup[i], std::vector<double>() ) );
}


allMasses::allMasses(unsigned int inMasses, unsigned int inEvents,
		std::vector<std::pair<unsigned int, unsigned int> >& inTup) :
				nInvMasses(inMasses),nEvents(inEvents){
	//alocate memory in advance
	for(unsigned int i=0; i<inTup.size(); i++)
		masses_sq.insert( std::make_pair( inTup[i], std::vector<double>(inEvents,0.) ) );
	eff = std::vector<double>(nEvents,1.);
	weight = std::vector<double>(nEvents,1.);
}
double allMasses::getReWeight(){ return reWeight; };

allMasses::allMasses():nInvMasses(0),nEvents(0) {}

void allMasses::resetWeights(){
	weight = std::vector<double>(nEvents,1.);
}
void allMasses::setEfficiency(double constEff){
	unsigned int nEvents = masses_sq.at(std::make_pair(2,3)).size();
	eff = std::vector<double>(nEvents,constEff);
}
void allMasses::setEfficiency(std::shared_ptr<Efficiency> effObj){
	unsigned int nEvents = masses_sq.at(std::make_pair(2,3)).size();
	eff = std::vector<double>(nEvents,1.);
	for(unsigned int i=0; i<nEvents;i++){
		std::vector<double> data;
		data.push_back(masses_sq.at( std::make_pair(2,3) )[i]);
		data.push_back(masses_sq.at( std::make_pair(1,3) )[i]);
		double value  = effObj->evaluate(data);
		/* We need to use sqrt(eff) here because in the current
		 * implementation the Amplitude value is squared after
		 * multiplication with the efficiency */
		eff.at(i) = sqrt(value);
		if(value==0) eff.at(i) = 0.001;
	}
	return;
}
