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

Kinematics* Kinematics::instance()
{
	if(!_inst) {
		throw std::runtime_error("No instance of Kinematics created! "
				"Create one first!");
	}

	return Kinematics::_inst;
}

Kinematics* Kinematics::_inst = 0;

double Kinematics::qSqValue(double sqrtS, double ma, double mb)
{
	double mapb = ma + mb;
	double mamb = ma - mb;
	double xSq = sqrtS*sqrtS;
	double t1 = xSq - mapb*mapb;
	double t2 = xSq - mamb*mamb;
	return ( t1*t2/(4*xSq) );
}

std::complex<double> Kinematics::qValue(double sqrtS, double ma, double mb)
{
	return Kinematics::phspFactor(sqrtS,ma,mb)*8.0*M_PI*sqrtS;
}

double Kinematics::FormFactor(double sqrtS, double ma, double mb, double spin,
		double mesonRadius,	formFactorType type)
{
	if(type == formFactorType::BlattWeisskopf && spin == 0){ return 1.0; }

	std::complex<double> qValue = Kinematics::qValue(sqrtS,ma,mb);
	return Kinematics::FormFactor(sqrtS, ma, mb, spin, mesonRadius, qValue, type);
}

double Kinematics::FormFactor(double sqrtS, double ma, double mb, double spin,
		double mesonRadius,	std::complex<double> qValue, formFactorType type)
{
	if(mesonRadius==0) return 1; //disable form factors
	if(type == formFactorType::noFormFactor) return 1; //disable form factors
	if(type == formFactorType::BlattWeisskopf && spin == 0){ return 1.0; }

	//From factor for a0(980) used by Crystal Barrel Phys.Rev.D78-074023
	if(type == formFactorType::CrystalBarrel){
		if (spin == 0) {
			double qSq = std::norm(qValue);
			double alpha = mesonRadius*mesonRadius/6;
			return std::exp(-alpha*qSq);
		}
		else
			throw std::runtime_error("Kinematics::FormFactor() | "
					"Form factors of type "
					+std::string(formFactorTypeString[type])
		+ " are implemented for spin 0 only!"
			);
	}

	//Blatt-Weisskopt form factors with normalization F(x=mR) = 1.
	//Reference: S.U.Chung Annalen der Physik 4(1995) 404-430
	//z = q / (interaction range). For the interaction range we assume 1/mesonRadius
	if(type == formFactorType::BlattWeisskopf){
		if (spin == 0) return 1;
		double qSq = std::norm(qValue);
		double z = qSq*mesonRadius*mesonRadius;
		/* Events below threshold
		 * What should we do if event is below threshold? Shouldn't really influence the result
		 * because resonances at threshold don't have spin(?) */
		z = std::fabs(z);

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
			throw std::runtime_error("Kinematics::FormFactor() | Form factors of type "
					+std::string(formFactorTypeString[type]) + " are implemented for spins up to 4!");
	}
	throw std::runtime_error("Kinematics::FormFactor() | Form factor type "
			+std::to_string((long long int)type)+" not specified!");

	return 0;
}

std::complex<double> Kinematics::phspFactor(double sqrtS, double ma, double mb)
{
	double s = sqrtS*sqrtS;
	std::complex<double> i(0,1);

	// == Two types of analytic continuation
	// 1) Complex sqrt
	//	std::complex<double> rhoOld;
	//	rhoOld = sqrt(std::complex<double>(qSqValue(sqrtS,ma,mb))) / (8*M_PI*sqrtS); //PDG definition
	//	rhoOld = sqrt(std::complex<double>(qSqValue(sqrtS,ma,mb))) / (0.5*sqrtS); //BaBar definition
	//	return rhoOld; //complex sqrt

	/* 2) Correct analytic continuation
	 * proper analytic continuation (PDG 2014 - Resonances (47.2.2))
	 * I'm not sure of this is correct for the case of two different masses ma and mb.
	 * Furthermore we divide by the factor 16*Pi*Sqrt[s]). This is more or less arbitrary
	 * and not mentioned in the reference, but it leads to a good agreement between both
	 * approaches.
	 */
	std::complex<double> rho;
	double q = std::sqrt( std::fabs(Kinematics::qSqValue(sqrtS,ma,mb)*4/s) );
	if( s<=0 ){ //below 0
		rho = -q/M_PI*std::log(std::fabs((1+q)/(1-q)));
	} else if( 0<s && s<=(ma+mb)*(ma+mb) ){ //below threshold
		rho = ( -2.0*q/M_PI * atan(1/q) ) / (i*16.0*M_PI*sqrtS);
	} else if( (ma+mb)*(ma+mb)<s ){//above threshold
		rho = ( -q/M_PI*log( std::fabs((1+q)/(1-q)) )+i*q ) / (i*16.0*M_PI*sqrtS);
	} else
		throw std::runtime_error("Kinematics::phspFactor| PhspFactor not "
				"defined at sqrtS = "+std::to_string((long double)sqrtS));

	if(rho.real()!=rho.real() || rho.imag()!=rho.imag())
		throw std::runtime_error("Kinematics::phspFactor| Result invalid (NaN)!");

	return rho; //correct analytical continuation
}

unsigned int Kinematics::FindVariable(std::string varName) const
{
	for(unsigned int i=0; i< _varNames.size(); i++)
		if(_varNames.at(i)==varName) return i;
	throw BadParameter("Kinematics::findVariable() | "
			"Variable"+varName+" not found!");
	return -999;
}

//======================== TWO-BODY KINEMATICS =============================
TwoBodyKinematics::TwoBodyKinematics(std::string _nameMother,
		std::string _name1, std::string _name2,	double deltaMassWindow) :
				name1(_name1), name2(_name2), Kinematics(_nameMother,0.0,2)
{
	_M = PhysConst::instance()->getMass(_nameMother);
	m1 = PhysConst::instance()->getMass(_name1);
	m2 = PhysConst::instance()->getMass(_name2);
	_spinM = PhysConst::instance()->getJ(_nameMother);
	spin1 = PhysConst::instance()->getJ(_name1);
	spin2 = PhysConst::instance()->getJ(_name2);
	if(_M==-999 || m1==-999|| m2==-999)
		throw std::runtime_error("TwoBodyKinematics(): Masses not set!");
	mass_min=((_M-deltaMassWindow)); mass_max=((_M+deltaMassWindow));
	mass_sq_max = mass_max*mass_max;
	mass_sq_min = mass_min*mass_max;
	_varNames.push_back("msq");

	init();
}

void TwoBodyKinematics::init()
{

}

bool TwoBodyKinematics::IsWithinPhsp(const dataPoint& point) const
{
	return 1;
	if(point.getVal(0)>=mass_sq_min && point.getVal(0)<=mass_sq_max) return 1;
	return 0;
}

void TwoBodyKinematics::EventToDataPoint(const Event& ev, dataPoint& point) const
{
	double weight = ev.getWeight();
	point.setWeight(weight);//reset weight
	const Particle& part1 = ev.getParticle(0);
	const Particle& part2 = ev.getParticle(1);
	double msq = Particle::invariantMass(part1,part2);
	point.setVal(0,msq);
	return;
}

//! get mass of particles
double TwoBodyKinematics::GetMass(unsigned int num) const
{
	if(num==0) return _M;
	if(num==1) return m1;
	if(num==2) return m2;
	throw std::runtime_error("TwoBodyKinematics::getMass(int) | "
			"Wrong particle requested!");
	return -999;
}

//! get mass of paticles
double TwoBodyKinematics::GetMass(std::string name) const
{
	if(name==_nameMother) return _M;
	if(name==name1) return m1;
	if(name==name2) return m2;
	throw std::runtime_error("TwoBodyKinematics::getMass(int) | "
			"Wrong particle "+name+" requested!");
	return -999;
}
