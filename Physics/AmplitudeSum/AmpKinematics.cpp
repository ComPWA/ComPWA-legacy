//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

#include "Physics/AmplitudeSum/AmpKinematics.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

AmpKinematics::AmpKinematics(std::shared_ptr<DoubleParameter> mR, int subSys, int spin, int m, int n,
		std::shared_ptr<DoubleParameter> mesonRadius, std::shared_ptr<DoubleParameter> motherRadius) :
		_M(-999), _mR(mR),_subSys(subSys), _spin(spin),_m(m),_n(n),
		_mesonRadius(mesonRadius), _motherRadius(motherRadius), _wignerD(subSys,spin)
{
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	_M=kin->M;
	if(_subSys==5){
		_ma=kin->m3;
		_mb=kin->m2;
		_mc=kin->m1;}
	if(_subSys==4){
		_ma=kin->m3;
		_mb=kin->m1;
		_mc=kin->m2;}
	if(_subSys==3){
		_ma=kin->m2;
		_mb=kin->m1;
		_mc=kin->m3;}
}

double AmpKinematics::qSqValue(double sqrtS, double ma, double mb){
	double mapb = ma + mb;
	double mamb = ma - mb;
	double xSq = sqrtS*sqrtS;
	double t1 = xSq - mapb*mapb;
	double t2 = xSq - mamb*mamb;
	return ( t1*t2/(4*xSq) );
}
std::complex<double> AmpKinematics::qValue(double sqrtS, double ma, double mb){
	//std::complex<double> result( sqrt( qSqValue(sqrtS,ma,mb) ) ); //complex sqrt!
	std::complex<double> result = phspFactor(sqrtS,ma,mb)*8.0*M_PI*sqrtS; //calculate from phsp factor
	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"AmpKinematics::qValue() | NaN! sqrtS="<<sqrtS<<" ma="<<ma<<" mb="<<mb<<std::endl;
	}
	return result;
}

double AmpKinematics::FormFactor(double sqrtS,double ma, double mb, double spin, double mesonRadius){
	if (spin == 0) return 1;
	//Blatt-Weisskopt form factors with normalization F(x=mR) = 1.
	//Reference: S.U.Chung Annalen der Physik 4(1995) 404-430
	//z = q / (interaction range). For the interaction range we assume 1/mesonRadius
	double z = qSqValue(sqrtS,ma,mb)*mesonRadius*mesonRadius;
	/* Events below threshold
	 * What should we do if event is below threshold? Shouldn't influence in
	 * practise because resonances at threshold don't have spin(?) */
	z = abs(z);

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
	else{
		std::cout<<"Wrong spin value! BLW factors only implemented for spin 0-4! "<<std::endl;
	}
	return 0;
}

std::complex<double> AmpKinematics::phspFactor(double sqrtS, double ma, double mb){
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
	double q = std::sqrt( std::abs(qSqValue(sqrtS,ma,mb)*4/s) );
	if( s<=0 ){ //below 0
		rho = -q/M_PI*std::log(std::abs((1+q)/(1-q)));
	} else if( 0<s && s<=(ma+mb)*(ma+mb) ){ //below threshold
		rho = ( -2.0*q/M_PI * atan(1/q) ) / (i*16.0*M_PI*sqrtS);
	} else if( (ma+mb)*(ma+mb)<s ){//above threshold
		rho = ( -q/M_PI*log( abs((1+q)/(1-q)) )+i*q ) / (i*16.0*M_PI*sqrtS);
	} else
		throw std::runtime_error("AmpKinematics::phspFactor| phspFactor not defined at sqrtS = "
				+std::to_string((long double)sqrtS));

	if(rho.real()!=rho.real() || rho.imag()!=rho.imag()){
		std::cout<<"AmpKinematics::phspFactor() | NaN! sqrtS="<<sqrtS<<" ma="<<ma<<" mb="<<mb<<std::endl;
	}
	return rho; //correct analytical continuation
}

std::complex<double> AmpKinematics::widthToCoupling(double mSq, double mR, double width,
		double ma, double mb, double spin, double mesonRadius)
{
	double sqrtS = sqrt(mSq);
	//calculate gammaA(s_R)
	double ffR = FormFactor(mR,ma,mb,spin,mesonRadius);
	std::complex<double> qR = qValue(mR,ma,mb);
	//calculate phsp factor
	std::complex<double> rho = phspFactor(sqrtS,ma,mb);
	std::complex<double> denom = std::pow(qR,spin)*ffR*sqrt(rho);
	std::complex<double> result = std::complex<double>(sqrt(mR*width), 0) / denom;
	return result;
}

std::complex<double> AmpKinematics::couplingToWidth(double mSq, double mR, double g,
		double ma, double mb, double spin, double mesonRadius)
{
	double sqrtM = sqrt(mSq);
	//calculate gammaA(s_R)
	double ffR = FormFactor(mR,ma,mb,spin,mesonRadius);
	std::complex<double> qR = std::pow(qValue(mR,ma,mb),spin);
	std::complex<double> gammaA = ffR*qR;
	//calculate phsp factor
	std::complex<double> rho = phspFactor(sqrtM,ma,mb);
	std::complex<double> result = std::norm(gammaA)*g*g*rho/ mR;
	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"AmpKinematics::couplingToWidth() | NaN! mSq="<<mSq
				<<" mR="<<mR<<" g="<<g<<" ma="<<ma<<" mb="<<mb<<std::endl;
		std::cout<<qR<<" "<<gammaA<<" "<<rho<<" "<<g<<std::endl;
	}
	return result;
}
