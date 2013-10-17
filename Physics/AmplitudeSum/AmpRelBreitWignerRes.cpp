//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "RooRealVar.h"

AmpRelBreitWignerRes::AmpRelBreitWignerRes(const char *name, const char *title,
		RooAbsReal& resMass, RooAbsReal& resWidth,
		RooAbsReal& d, //  meson radius
		int subSys,
		Int_t resSpin,
		Int_t m,
		Int_t n
		) :
		AmpAbsDynamicalFunction(name,title),
		AmpKinematics(resMass.getVal(), subSys, resSpin, m, n, AmpKinematics::barrierType(BWPrime), d.getVal(), 1.5),
//		_x13("x13", "Observable", this, 0),
//		_x23("x23", "Observable",this, 0),
		_resWidth("resWidth", "Width", this, resWidth),
//		_subSys(subSys),
		_wignerD(name, title, (UInt_t)resSpin, (UInt_t)m, (UInt_t)n, subSys)
{
	initialise();
}
AmpRelBreitWignerRes::AmpRelBreitWignerRes(const char *name, const char *title,
		RooAbsReal& x13, //  mass at which to evaluate RBW
		RooAbsReal& x23, //  mass at which to evaluate RBW
		RooAbsReal& resMass, RooAbsReal& resWidth,
		RooAbsReal& d, //  meson radius
		int subSys,
		Int_t resSpin,
		Int_t m,
		Int_t n
		) :
		AmpAbsDynamicalFunction(name,title),
		AmpKinematics(resMass.getVal(), subSys, resSpin, m, n, AmpKinematics::barrierType(BWPrime), d.getVal(), 1.5),
		_x13("x13", "Observable", this, x13),
		_x23("x23", "Observable",this, x23),
		_resWidth("resWidth", "Width", this, resWidth),
//		_subSys(subSys),
		_wignerD(name, title, x13, x23, (UInt_t)resSpin, (UInt_t)m, (UInt_t)n, subSys)
{
	initialise();
}


AmpRelBreitWignerRes::AmpRelBreitWignerRes(const AmpRelBreitWignerRes& other, const char* newname) :
		  AmpAbsDynamicalFunction(other, newname),
		  AmpKinematics(other),
		  _x13("x13", this, other._x13),
		  _x23("x23", this, other._x23),
		  _resWidth("resWidth", this, other._resWidth),
		  _wignerD(other._wignerD)
{
	initialise();
}

AmpRelBreitWignerRes::AmpRelBreitWignerRes(const AmpRelBreitWignerRes& other) :
		  AmpAbsDynamicalFunction(other.GetName(), other.GetTitle()),
		  AmpKinematics(other),
		  _x13("x13", this, other._x13),
		  _x23("x23", this, other._x23),
		  _resWidth("resWidth", this, other._resWidth),
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
RooComplex AmpRelBreitWignerRes::evaluate() const {
	if(_ma == -999 || _mb == -999 ||_mc == -999 ||_M == -999){
		std::cout<<"Masses of decay products not set!"<<std::endl;
		return 0;
	}
	RooComplex result;
//	double m  = Double_t(_x23);
	double m = dataPoint::instance()->getM(_subSys);
	double spinTerm = _wignerD.evaluate();
	double Gamma0 = double(_resWidth);
	double GammaV = Gamma0 * pow(q(m) / q0(), 2.*_spin + 1.) * (_mR / m) * BLres2(m);

	RooComplex denom = RooComplex(_mR*_mR - m*m, -_mR * GammaV);

//	  result = RooComplex(spinTerm*_mR * Gamma0) / denom; //wrong!
	//  result = RooComplex(spinTerm*BLprime2(m)) / denom;
	result = RooComplex( spinTerm ) / denom; //Laura++ (old) definition - is this used in DKsKK analysis?
//	result = RooComplex( spinTerm*sqrt(BLres2(m))*sqrt(BLmother2(m)) ) / denom; //Laura++ (new) definition

	if(result.re()!=result.re()) {std::cout << "RE part NAN" << std::endl;return 0;}
	if(result.im()!=result.im()) {std::cout << "IM part NAN" << std::endl; return 0;}
	return result;
}
//double AmpRelBreitWignerRes::evaluate(double x[], int dim, void * param) const {
//
//	if( !_dpPoint.DPKin.isWithinDP(x[0],x[1]) ) return 0;
//
//	return 1;
//}

////  implement stubs for virtual functions
////  to satisfy interface of RooAbsArg.
////  not sure which ones are really needed.

TObject*  AmpRelBreitWignerRes::clone (const char *newname) const 
{
	std::cout << "Unimplemented AmpRelBreitWignerRes::clone called for " << this << std::endl
			<< "This will crash the program! "
			<< std::endl;

	//  cout << __PRETTY_FUNCTION__ << "(" << ((newname)?newname:"(null)") << ")" << endl;
	return 0;//new AmpRelBreitWignerRes(*this, newname);
}

Bool_t AmpRelBreitWignerRes::readFromStream(std::istream&, Bool_t, Bool_t) {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	return false;
}

void AmpRelBreitWignerRes::writeToStream(std::ostream&, Bool_t) const {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	;
}

Bool_t AmpRelBreitWignerRes::operator==(const RooAbsArg&) {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	return false;
}

void AmpRelBreitWignerRes::syncCache(const RooArgSet*)  {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	;
}

void AmpRelBreitWignerRes::copyCache(const RooAbsArg*, Bool_t, Bool_t) {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	;
}

void AmpRelBreitWignerRes::attachToTree(TTree&, Int_t) {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	;
}

void AmpRelBreitWignerRes::attachToVStore(RooVectorDataStore&){
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	;
}

void AmpRelBreitWignerRes::setTreeBranchStatus(TTree&, Bool_t) {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	;
}

void AmpRelBreitWignerRes::fillTreeBranch(TTree&) {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	;
}

bool AmpRelBreitWignerRes::isIdentical(const RooAbsArg&, Bool_t){
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	return 0;
}

RooAbsArg* AmpRelBreitWignerRes::createFundamental(const char*) const {
	std::cout << "Unimplemented " << __PRETTY_FUNCTION__ << " called for " << this << std::endl
			<< "This might crash the program! "
			<< std::endl;
	return 0;
}
