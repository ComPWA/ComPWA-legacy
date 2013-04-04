//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "RooRealVar.h"

AmpGausRes::AmpGausRes(const char *name, const char *title,
					   RooAbsReal& x, ///  mass at which to evaluate RBW
					   RooAbsReal& resMass, RooAbsReal& resWidth,
					   int subSys) :
  AmpAbsDynamicalFunction(name,title),
  _x("x", "Observable", this, x),
  _m0("m0", "Mass", this, resMass),
  _resWidth("resWidth", "Width", this, resWidth),
  _subSys(subSys)
{
  initialise();
}


AmpGausRes::AmpGausRes(const AmpGausRes& other, const char* newname) :
  AmpAbsDynamicalFunction(other, newname),
  _x("x", this, other._x),
  _m0("m0", this, other._m0),
  _resWidth("resWidth", this, other._resWidth),
  _subSys(other._subSys)
{
  initialise();
}

AmpGausRes::AmpGausRes(const AmpGausRes& other) :
  AmpAbsDynamicalFunction(other.GetName(), other.GetTitle()),
  _x("x", this, other._x),
  _m0("m0", this, other._m0),
  _resWidth("resWidth", this, other._resWidth),
  _subSys(other._subSys)
{
  initialise();
}

AmpGausRes::~AmpGausRes() 
{
}

void AmpGausRes::initialise() 
{
}   

RooComplex AmpGausRes::evaluate() const {


  double m0 = Double_t(_m0);
  double m  = Double_t(_x);
  
  RooComplex gaus = RooComplex(exp(-1*(m-m0)*(m-m0)/_resWidth/_resWidth/2.),0);

  return gaus;
}
                       
////  implement stubs for virtual functions
////  to satisfy interface of RooAbsArg.
////  not sure which ones are really needed.

TObject*  AmpGausRes::clone (const char *newname) const 
{
  //  std::cout << "Unimplemented AmpRelBreitWignerRes::clone called for " << this << endl
  //	    << "This will crash the program! " 
  //	    << endl;

  //  cout << __PRETTY_FUNCTION__ << "(" << ((newname)?newname:"(null)") << ")" << endl;
  return new AmpGausRes(*this, newname);
}

Bool_t AmpGausRes::readFromStream(std::istream&, Bool_t, Bool_t) {
  cout << __PRETTY_FUNCTION__ << endl;
  return false;
}

void AmpGausRes::writeToStream(std::ostream&, Bool_t) const {
  cout << __PRETTY_FUNCTION__ << endl;
 ; 
}

Bool_t AmpGausRes::operator==(const RooAbsArg&) {
  cout << __PRETTY_FUNCTION__ << endl;
  return false;
}

void AmpGausRes::syncCache(const RooArgSet*)  {
  cout << __PRETTY_FUNCTION__ << endl;
  ;
}

void AmpGausRes::copyCache(const RooAbsArg*, Bool_t, Bool_t) {
  cout << __PRETTY_FUNCTION__ << endl;
  ; 
}

void AmpGausRes::attachToTree(TTree&, Int_t) {
  cout << __PRETTY_FUNCTION__ << endl;
  ;
}

void AmpGausRes::attachToVStore(RooVectorDataStore&) {
  cout << __PRETTY_FUNCTION__ << endl;
  ;
}

void AmpGausRes::setTreeBranchStatus(TTree&, Bool_t) {
  cout << __PRETTY_FUNCTION__ << endl;
  ;
}

void AmpGausRes::fillTreeBranch(TTree&) {
  cout << __PRETTY_FUNCTION__ << endl;
  ;
}
RooAbsArg* AmpGausRes::createFundamental(const char*) const {
  std::cout << "Unimplemented " << __PRETTY_FUNCTION__ << " called for " << this << endl
	    << "This might crash the program! " 
	    << endl;
  return 0;
}
