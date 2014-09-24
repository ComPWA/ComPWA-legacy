//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include <stdlib.h>


AmpRelBreitWignerRes::AmpRelBreitWignerRes(const char *name,
		DoubleParameter& resMass, DoubleParameter& resWidth,
		DoubleParameter& mesonRadius, DoubleParameter& motherRadius,
		int subSys,
		int resSpin, int m, int n) :
		AmpAbsDynamicalFunction(name),
		AmpKinematics(resMass, subSys, resSpin, m, n, AmpKinematics::barrierType(BWPrime),
				mesonRadius, motherRadius),
		_resWidth(resWidth),
		//_wignerD(name, resSpin, m, n, subSys)
		_wignerD(subSys, resSpin),
		foundMasses(false),
		nParams(6)
{
	initialise();
}

AmpRelBreitWignerRes::AmpRelBreitWignerRes(const AmpRelBreitWignerRes& other, const char* newname) :
														  AmpAbsDynamicalFunction(other, newname),
														  AmpKinematics(other),
														  _resWidth(other._resWidth),
														  _wignerD(other._wignerD)
{
	initialise();
}

AmpRelBreitWignerRes::AmpRelBreitWignerRes(const AmpRelBreitWignerRes& other) :
														  AmpAbsDynamicalFunction(other),
														  AmpKinematics(other),
														  _resWidth(other._resWidth),
														  _wignerD(other._wignerD)
{
	initialise();
}

AmpRelBreitWignerRes::~AmpRelBreitWignerRes() 
{
}

void AmpRelBreitWignerRes::initialise() 
{
}

void AmpRelBreitWignerRes::setDecayMasses(double ma, double mb, double mc, double M){
	AmpKinematics::setDecayMasses(ma, mb, mc, M);
	//	_wignerD.setDecayMasses(ma, mb, mc, M);
	return;
}
std::complex<double> AmpRelBreitWignerRes::evaluateAmp(dataPoint& point) {
	if(_ma == -999 || _mb == -999 ||_mc == -999 ||_M == -999){
		std::cout<<"Masses of decay products not set!"<<std::endl;
		return 0;
	}
	std::complex<double> result;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	if(!foundMasses){
		id23 = point.getID("m23sq");
		id13 = point.getID("m13sq");
		foundMasses=true;
	}
	double m23sq = point.getVal(id23);
	double m13sq = point.getVal(id13);
	double m12sq = kin->getThirdVariableSq(m23sq,m13sq);
	double mSq = -999;
	switch(_subSys){
	case 3: mSq=(m12sq); break;
	case 4: mSq=(m13sq); break;
	case 5: mSq=(m23sq); break;
	}
	//		double BLWeiss2 = BLres2(sqrt(mSq));
	//		double qTerm = std::pow((q(sqrt(mSq)) / q0()), (2.*_spin + 1.));
	//		double Gamma0 = _resWidth.GetValue();
	//		double GammaV = Gamma0 * qTerm * (_mR / sqrt(mSq)) * BLWeiss2;
	//		std::complex<double> denom(_mR*_mR - mSq, -_mR * GammaV);
	//
	//		if(sqrt(m23sq)==2.84515) std::cout << " DEBUG2 " << q(m) << " " << q0() << std::endl;
	//
	//		result = std::complex<double>( _norm ) / denom; //Laura++ (old) definition
	//	//	result = std::complex<double>(_norm*sqrt(BLWeiss2)) / denom;
	//	//	result = std::complex<double>( _norm*sqrt(BLWeiss2)*sqrt(BLmother2(m)) ) / denom; //Laura++ (new) definition
	//
	//		if(result.real()!=result.real()) {std::cout << "RE part NAN" << std::endl;return 0;}
	//		if(result.imag()!=result.imag()) {std::cout << "IM part NAN" << std::endl; return 0;}
	//if(_mR==0.783) {std::cout << "Omega Norm " << _norm << std::endl; return 0;}
	//		return result;

	//	return (dynamicalFunction(mSq,_mR,_ma,_mb,_resWidth.GetValue(),_spin,_mesonRadius)*_norm);
	return dynamicalFunction(mSq,_mR.GetValue(),_ma,_mb,_resWidth.GetValue(),_spin,_mesonRadius);
}
std::complex<double> AmpRelBreitWignerRes::dynamicalFunction(double mSq, double mR, double ma, double mb, double gamma0, unsigned int J, double mesonRadius){
	double m = sqrt(mSq);
	double BLWeiss2 = AmpKinematics::BlattWeiss(m,mR,ma,mb,J,mesonRadius);
	double qTerm = std::pow((qValue(m,ma,mb) / qValue(mR,ma,mb)), (2.*J+ 1.));
	double GammaV = gamma0 * qTerm * (mR / m) * BLWeiss2;
	std::complex<double> denom(mR*mR - m*m, -mR * GammaV);

	std::complex<double> result = std::complex<double>( 1. ,0 ) / denom; //Laura++ (old) definition
	//	result = std::complex<double>( (2.*J+1.)*sqrt(BLWeiss2) ) / denom;
	//	result = std::complex<double>( (2.*J+1.)*sqrt(BLWeiss2)*sqrt(BLmother2(m)) ) / denom; //Laura++ (new) definition

	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"nan in BW: "<<BLWeiss2<<" "<<mR<<" "<<mSq<<" "<<ma<<" "<<mb<<std::endl;
		return 0;
	}
	return result;
}

/*std::complex<double> AmpRelBreitWignerRes::evaluateTree(const ParameterList& paras) const {
    double Gamma0, GammaV;
    double m0 = double(paras.GetParameterValue("m0_"+_name));
    double m  = double(paras.GetParameterValue("x_"+_name));
    double ma = double(paras.GetParameterValue("ma_"+_name));
    double mb = double(paras.GetParameterValue("mb_"+_name));
    unsigned int spin = double(paras.GetParameterValue("spin_"+_name));
    double d = double(paras.GetParameterValue("d_"+_name));
    double norm = double(paras.GetParameterValue("norm_"+_name));

    Gamma0 = double(paras.GetParameterValue("resWidth_"+_name));
   // GammaV = Gamma0 * (m0 / m) * pow(q(ma,mb,m) / q0(ma,mb,m0), 2.*spin + 1.)  * BLprime2(ma,mb,m0,m,d,spin);
	double BLWeiss2 = BLres2(m);
	double qTerm = pow(q(m) / q0(), 2.*spin + 1.);
	//double Gamma0 = _resWidth.GetValue();
	GammaV = Gamma0 * qTerm * (m0 / m) * BLWeiss2;
	std::complex<double> denom(m0*m0 - m*m, -m0 * GammaV);

    //RooComplex denom = RooComplex(m0*m0 - m*m, -m0 * GammaV);
    //RooComplex res = RooComplex(m0 * Gamma0) / denom; //TODO: use same functions as over evaluates

    std::complex<double> result = std::complex<double>( norm ) / denom;
    //std::shared_ptr<ComplexParameter> bw(new ComplexParameter("relBW of "+name, result));
    return result*_wignerD.evaluateTree(paras,_name);
}*/
