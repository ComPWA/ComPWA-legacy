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

#include <boost/log/trivial.hpp>
using namespace boost::log;

#include "Core/DataPoint.hpp"
#include "Core/Kinematics.hpp"
#include "Core/PhysConst.hpp"

Kinematics* Kinematics::instance(){
	if(!_inst) {
		throw std::runtime_error("No instance of Kinematics created! Create one first!");
	}

	return Kinematics::_inst;
}

Kinematics* Kinematics::_inst = 0;

TwoBodyKinematics::TwoBodyKinematics(std::string _nameMother, std::string _name1, std::string _name2,
		double deltaMassWindow) {
	nPart = 2;
	M = PhysConst::instance()->getMass(_nameMother);
	m1 = PhysConst::instance()->getMass(_name1);
	m2 = PhysConst::instance()->getMass(_name2);
	spinM = PhysConst::instance()->getJ(_nameMother);
	spin1 = PhysConst::instance()->getJ(_name1);
	spin2 = PhysConst::instance()->getJ(_name2);
	if(M==-999 || m1==-999|| m2==-999)
		throw std::runtime_error("TwoBodyKinematics(): Masses not set!");
	mass_min=((M-deltaMassWindow)); mass_max=((M+deltaMassWindow));
	mass_sq_max = mass_max*mass_max;
	mass_sq_min = mass_min*mass_max;
	varNames.push_back("msq");

	init();
}
void TwoBodyKinematics::init(){
}
bool TwoBodyKinematics::isWithinPhsp(const dataPoint& point){
	return 1;
	if(point.getVal(0)>=mass_sq_min && point.getVal(0)<=mass_sq_max) return 1;
	return 0;
}
void TwoBodyKinematics::eventToDataPoint(Event& ev, dataPoint& point){
	double weight = ev.getWeight();
	point.setWeight(weight);//reset weight
	Particle part1 = ev.getParticle(0);
	Particle part2 = ev.getParticle(1);
	double msq = Particle::invariantMass(part1,part2);
	point.setVal(0,msq);
	return;
}
//! get mass of particles
double TwoBodyKinematics::getMass(unsigned int num){
	if(num==0) return M;
	if(num==1) return m1;
	if(num==2) return m2;
	throw std::runtime_error("TwoBodyKinematics::getMass(int) | wrong particle requested!");
	return -999;
}
//! get mass of paticles
double TwoBodyKinematics::getMass(std::string name){
	if(name==nameMother) return M;
	if(name==name1) return m1;
	if(name==name2) return m2;
	throw std::runtime_error("TwoBodyKinematics::getMass(int) | wrong particle "+name+" requested!");
	return -999;
}

