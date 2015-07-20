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

double Kinematics::qSqValue(double sqrtS, double ma, double mb){
	double mapb = ma + mb;
	double mamb = ma - mb;
	double xSq = sqrtS*sqrtS;
	double t1 = xSq - mapb*mapb;
	double t2 = xSq - mamb*mamb;
	return ( t1*t2/(4*xSq) );
}
std::complex<double> Kinematics::qValue(double sqrtS, double ma, double mb){
	//std::complex<double> result( sqrt( qSqValue(sqrtS,ma,mb) ) ); //complex sqrt!
	std::complex<double> result = Kinematics::phspFactor(sqrtS,ma,mb)*8.0*M_PI*sqrtS; //calculate from phsp factor
	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"AmpKinematics::qValue() | NaN! sqrtS="<<sqrtS<<" ma="<<ma<<" mb="<<mb<<std::endl;
	}
	return result;
}

double Kinematics::FormFactor(double sqrtS, double ma, double mb, double spin, double mesonRadius){
	if (spin == 0) return 1;
	//Blatt-Weisskopt form factors with normalization F(x=mR) = 1.
	//Reference: S.U.Chung Annalen der Physik 4(1995) 404-430
	//z = q / (interaction range). For the interaction range we assume 1/mesonRadius
	double z = Kinematics::qSqValue(sqrtS,ma,mb)*mesonRadius*mesonRadius;
	/* Events below threshold
	 * What should we do if event is below threshold? Shouldn't really influence the result
	 * because resonances at threshold don't have spin(?) */
	z = std::abs(z);

	if (spin == 1){
		return( sqrt(2*z/(z+1)) );
	}
	else if (spin == 2) {
		return ( sqrt( 13*z*z/( (z-3)*(z-3)+9*z ) ) );
	}
	else if (spin == 3) {
		return ( sqrt( 277*z*z*z/( z*(z-15)*(z-15) + 9*(2*z-5) ) ) );
	}
	else if (spin == 4) {
		return ( sqrt( 12746*z*z*z*z/( (z*z-45*z+105)*(z*z-45*z+105) + 25*z*(2*z-21)*(2*z-21) ) ) );
	}
	else
		throw std::runtime_error("Kinematics::FormFactor() | Form factors only implemented for spins uo to 4!");
	std::cout<<"we should never reach this point!"<<std::endl;
	return 0;
}

std::complex<double> Kinematics::phspFactor(double sqrtS, double ma, double mb){
	double s = sqrtS*sqrtS;
	std::complex<double> i(0,1);
	std::complex<double> rho;
	std::complex<double> rhoOld;

	// == Two types of analytic continuation
	// 1) Complex sqrt
	//rhoOld = sqrt(std::complex<double>(qSqValue(sqrtS,ma,mb))) / (8*M_PI*sqrtS); //PDG definition
	//rhoOld = sqrt(std::complex<double>(qSqValue(sqrtS,ma,mb))) / (0.5*sqrtS); //BaBar definition
	//return rhoOld; //complex sqrt

	/* 2) Correct analytic continuation
	 * proper analytic continuation (PDG 2014 - Resonances (47.2.2))
	 * I'm not sure of this is correct for the case of two different masses ma and mb.
	 * Furthermore we divide by the factor 16*Pi*Sqrt[s]). This is more or less arbitrary
	 * and not mentioned in the reference, but it leads to a good agreement between both
	 * approaches.
	 */
	double q = std::sqrt( std::abs(Kinematics::qSqValue(sqrtS,ma,mb)*4/s) );
	if( s<=0 ){ //below 0
		rho = -q/M_PI*std::log(std::abs((1+q)/(1-q)));
	} else if( 0<s && s<=(ma+mb)*(ma+mb) ){ //below threshold
		rho = ( -2.0*q/M_PI * atan(1/q) ) / (i*16.0*M_PI*sqrtS);
	} else if( (ma+mb)*(ma+mb)<s ){//above threshold
		rho = ( -q/M_PI*log( std::abs((1+q)/(1-q)) )+i*q ) / (i*16.0*M_PI*sqrtS);
	} else
		throw std::runtime_error("AmpKinematics::phspFactor| phspFactor not defined at sqrtS = "
				+std::to_string((long double)sqrtS));

	if(rho.real()!=rho.real() || rho.imag()!=rho.imag()){
		std::cout<<"AmpKinematics::phspFactor() | NaN! sqrtS="<<sqrtS<<" ma="<<ma<<" mb="<<mb<<std::endl;
	}
	return rho; //correct analytical continuation
}
