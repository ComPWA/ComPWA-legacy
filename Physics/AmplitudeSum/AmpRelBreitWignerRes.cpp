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
		RooAbsReal& x1, //  mass at which to evaluate RBW
		RooAbsReal& resMass, RooAbsReal& resWidth,
		RooAbsReal& d, //  meson radius
		int subSys,
		Int_t resSpin,
		Int_t m,
		Int_t n,
		) :
		AmpAbsDynamicalFunction(name,title),
		AmpKinematics(resMass.getVal(), resSpin, AmpKinematics::barrierType(BWPrime), d.getVal(), 1.5),
		_x("x1", "Observable", this, x1),
		_resWidth("resWidth", "Width", this, resWidth),
		_subSys(subSys),
		_wignerD()
{
	initialise();
}


AmpRelBreitWignerRes::AmpRelBreitWignerRes(const AmpRelBreitWignerRes& other, const char* newname) :
		  AmpAbsDynamicalFunction(other, newname),
		  AmpKinematics(other),
		  _x("x", this, other._x),
		  _resWidth("resWidth", this, other._resWidth),
		  _subSys(other._subSys)
{
	initialise();
}

AmpRelBreitWignerRes::AmpRelBreitWignerRes(const AmpRelBreitWignerRes& other) :
		  AmpAbsDynamicalFunction(other.GetName(), other.GetTitle()),
		  AmpKinematics(other),
		  _x("x", this, other._x),
		  _resWidth("resWidth", this, other._resWidth),
		  _subSys(other._subSys)
{
	initialise();
}

AmpRelBreitWignerRes::~AmpRelBreitWignerRes() 
{
}

void AmpRelBreitWignerRes::initialise() 
{
}

RooComplex AmpRelBreitWignerRes::evaluate() const {
	if(_ma == -999 || _mb == -999 ||_mc == -999 ||_M == -999){
		std::cout<<"Masses of decay products not set!"<<std::endl;
		return 0;
	}
	RooComplex result;
	double m  = Double_t(_x);

	double Gamma0 = double(_resWidth);
	double GammaV = Gamma0 * pow(q(m) / q0(), 2.*_spin + 1.) * (_mR / m) * BLres2(m);
//	std::cout<<q(m,_ma,_mb) << " "<<q0(_ma,_mb)<<std::endl;

	RooComplex denom = RooComplex(_mR*_mR - m*m, -_mR * GammaV);

//	  result = RooComplex(_mR * Gamma0) / denom; //wrong!
	//  result = RooComplex(BLprime2(m)) / denom;
	result = RooComplex( 1 ) / denom; //Laura++ (old) definition - is this used in DKsKK analysis?
//	result = RooComplex( sqrt(BLres2(m))*sqrt(BLmother2(m)) ) / denom; //Laura++ (new) definition
//	  std::cout<<"-- "<<_resWidth<<" " <<Gamma0<<" "<<_mR<<std::endl;
//	  std::cout<<sqrt(BLres2(m))<< " "<<sqrt(BLmother2(m))<< " "<<_mR*Gamma0<<std::endl;

	if(result.re()!=result.re()) {std::cout << "RE part NAN" << std::endl; return 0;}
	if(result.im()!=result.im()) {std::cout << "IM part NAN" << std::endl; return 0;}
	return result;
}

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
