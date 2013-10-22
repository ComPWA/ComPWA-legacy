//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"

AmpRelBreitWignerRes::AmpRelBreitWignerRes(const char *name,
		DoubleParameter& resMass, DoubleParameter& resWidth,
		double& mesonRadius, //  meson radius
		int subSys,
		int resSpin, int m, int n
		) :
		AmpAbsDynamicalFunction(name),
		AmpKinematics(resMass, subSys, resSpin, m, n, AmpKinematics::barrierType(BWPrime), mesonRadius, 1.5),
		_resWidth(resWidth),
		_wignerD(name, resSpin, m, n, subSys)
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
	  _wignerD.setDecayMasses(ma, mb, mc, M);
	  return;
  }
std::complex<double> AmpRelBreitWignerRes::evaluate() const {
	if(_ma == -999 || _mb == -999 ||_mc == -999 ||_M == -999){
		std::cout<<"Masses of decay products not set!"<<std::endl;
		return 0;
	}
	std::complex<double> result;
	double m = dataPoint::instance()->getM(_subSys);
	double spinTerm = _wignerD.evaluate();
//	double mr = _mR.GetValue();
	double Gamma0 = _resWidth.GetValue();
	double GammaV = Gamma0 * pow(q(m) / q0(), 2.*_spin + 1.) * (_mR / m) * BLres2(m);

	std::complex<double> denom(_mR*_mR - m*m, -_mR * GammaV);

//	  result = RooComplex(spinTerm*_mR * Gamma0) / denom; //wrong!
	//  result = RooComplex(spinTerm*BLprime2(m)) / denom;
	result = std::complex<double>( spinTerm ) / denom; //Laura++ (old) definition - is this used in DKsKK analysis?
//	result = RooComplex( spinTerm*sqrt(BLres2(m))*sqrt(BLmother2(m)) ) / denom; //Laura++ (new) definition

	if(result.real()!=result.real()) {std::cout << "RE part NAN" << std::endl;return 0;}
	if(result.imag()!=result.imag()) {std::cout << "IM part NAN" << std::endl; return 0;}
	return result;
}
//double AmpRelBreitWignerRes::evaluate(double x[], int dim, void * param) const {
//
//	if( !_dpPoint.DPKin.isWithinDP(x[0],x[1]) ) return 0;
//
//	return 1;
//}
