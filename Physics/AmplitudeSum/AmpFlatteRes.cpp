//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "RooRealVar.h"

AmpFlatteRes::AmpFlatteRes(const char *name, const char *title,
		RooAbsReal& x13, ///  mass at which to evaluate RBW
		RooAbsReal& x23, ///  mass at which to evaluate RBW
		RooAbsReal& resMass, RooAbsReal& resWidth,
		RooAbsReal& d, ///  meson radius
		RooAbsReal& par1,
		RooAbsReal& par2,
		int subSys, ///  meson radius
		Int_t resSpin,
		Int_t m,
		Int_t n) :
		AmpAbsDynamicalFunction(name,title),
		AmpKinematics(resMass.getVal(), resSpin, m, n, AmpKinematics::barrierType(BWPrime), d.getVal(), 1.5),
		_x13("x13", "Observable", this, x13),
		_x23("x23", "Observable", this, x23),
		_par1("par1","par1",this,par1),
		_par2("par2","par2",this,par2),
		_subSys(subSys),
		_wignerD(name, title, x13, x23, (UInt_t) resSpin,(UInt_t)  m,(UInt_t)  n)
{
	initialise();
}


AmpFlatteRes::AmpFlatteRes(const AmpFlatteRes& other, const char* newname) :
		  AmpAbsDynamicalFunction(other, newname),
		  AmpKinematics(other),
		  _x13("x13", "Observable", this, other._x13),
		  _x23("x23", "Observable", this, other._x23),
		  _par1("par1","par1",this,other._par1),
		  _par2("par2","par2",this,other._par2),
		  _subSys(other._subSys),
		  _wignerD(other._wignerD)
{
	initialise();
}

AmpFlatteRes::AmpFlatteRes(const AmpFlatteRes& other) :
		  AmpAbsDynamicalFunction(other.GetName(), other.GetTitle()),
		  AmpKinematics(other),
		  _x13("x13", "Observable", this, other._x13),
		  _x23("x23", "Observable", this, other._x23),
		  _par1("par1","par1",this,other._par1),
		  _par2("par2","par2",this,other._par2),
		  _subSys(other._subSys),
		  _wignerD(other._wignerD)
{
	initialise();
}

AmpFlatteRes::~AmpFlatteRes() 
{
}

void AmpFlatteRes::initialise() 
{
}
  void AmpFlatteRes::setDecayMasses(double ma, double mb, double mc, double M){
	  AmpKinematics::setDecayMasses(ma, mb, mc, M);
	  _wignerD.setDecayMasses(ma, mb, mc, M);
	  return;
  }

void AmpFlatteRes::setBarrierMass(double mBarA, double mBarB) {
	_mBarA = mBarA;
	_mBarB = mBarB;
}

RooComplex AmpFlatteRes::evaluate() const {
	if(_mBarA<0||_mBarA>5||_mBarB<0||_mBarB>5) {
		cout<<"Barrier masses not set! Use setBarrierMass() first!"<<endl;
		return 0;
	}
	//	double m0 = Double_t(_m0);
	double m  = Double_t(_x23);

	double p1 = 2*q(m, _mBarA,_mBarB)/m;//break-up momenta hidden channel (e.g. a0->eta pi)
	double p2 = 2*q(m, _ma,_mb)/m;//break-up momenta decay channel (e.g. a0->KK)
	double g1 = _par1;//couppling a0->eta pi
	double g2 = _par2;//coupling a0->KK

	RooComplex denom = RooComplex(_mR*_mR - m*m, -p1*g1*g1-p2*g2*g2);

	//	RooComplex result = (RooComplex(g2*g2,0) / denom); //use KK decay channel here
	RooComplex result = (RooComplex(g2,0) / denom); //use KK decay channel here

	if(result.re()!=result.re()) {std::cout << "RE part NAN" << std::endl; return 0;}
	if(result.im()!=result.im()) {std::cout << "IM part NAN" << std::endl; return 0;}
	return result;

}
////  implement stubs for virtual functions
////  to satisfy interface of RooAbsArg.
////  not sure which ones are really needed.

TObject*  AmpFlatteRes::clone (const char *newname) const 
{
	//  std::cout << "Unimplemented AmpRelBreitWignerRes::clone called for " << this << endl
	//	    << "This will crash the program! "
	//	    << endl;

	//  cout << __PRETTY_FUNCTION__ << "(" << ((newname)?newname:"(null)") << ")" << endl;
	return new AmpFlatteRes(*this, newname);
}

Bool_t AmpFlatteRes::readFromStream(std::istream&, Bool_t, Bool_t) {
	cout << __PRETTY_FUNCTION__ << endl;
	return false;
}

void AmpFlatteRes::writeToStream(std::ostream&, Bool_t) const {
	cout << __PRETTY_FUNCTION__ << endl;
	;
}

Bool_t AmpFlatteRes::operator==(const RooAbsArg&) {
	cout << __PRETTY_FUNCTION__ << endl;
	return false;
}

void AmpFlatteRes::syncCache(const RooArgSet*)  {
	cout << __PRETTY_FUNCTION__ << endl;
	;
}

void AmpFlatteRes::copyCache(const RooAbsArg*, Bool_t, Bool_t) {
	cout << __PRETTY_FUNCTION__ << endl;
	;
}

void AmpFlatteRes::attachToTree(TTree&, Int_t) {
	cout << __PRETTY_FUNCTION__ << endl;
	;
}

void AmpFlatteRes::attachToVStore(RooVectorDataStore&){
	cout << __PRETTY_FUNCTION__ << endl;
	;
}

void AmpFlatteRes::setTreeBranchStatus(TTree&, Bool_t) {
	cout << __PRETTY_FUNCTION__ << endl;
	;
}

void AmpFlatteRes::fillTreeBranch(TTree&) {
	cout << __PRETTY_FUNCTION__ << endl;
	;
}

RooAbsArg* AmpFlatteRes::createFundamental(const char*) const {
	std::cout << "Unimplemented " << __PRETTY_FUNCTION__ << " called for " << this << endl
			<< "This might crash the program! "
			<< endl;
	return 0;
}

Bool_t AmpFlatteRes::isIdentical(const RooAbsArg&, Bool_t){
	cout << __PRETTY_FUNCTION__ << endl;
	return 0;
}

