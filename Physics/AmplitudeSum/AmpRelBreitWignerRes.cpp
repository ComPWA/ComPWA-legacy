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
					   RooAbsReal& x, ///  mass at which to evaluate RBW
					   RooAbsReal& resMass, RooAbsReal& resWidth,
					   RooAbsReal& d, ///  meson radius
					   int subSys, ///  meson radius
					   Int_t resSpin) :
  AmpAbsDynamicalFunction(name,title),
  _x("x", "Observable", this, x),
  _m0("m0", "Mass", this, resMass),
  _resWidth("resWidth", "Width", this, resWidth),
  _d("d", "d", this, d),
  _subSys(subSys),
  _spin(resSpin)
{
  initialise();
}


AmpRelBreitWignerRes::AmpRelBreitWignerRes(const AmpRelBreitWignerRes& other, const char* newname) :
  AmpAbsDynamicalFunction(other, newname),
  _x("x", this, other._x),
  _m0("m0", this, other._m0),
  _resWidth("resWidth", this, other._resWidth),
  _d("d", this, other._d),
  _subSys(other._subSys),
  _spin(other._spin),
  _ma(other._ma),
  _mb(other._mb)
{
  initialise();
}

AmpRelBreitWignerRes::~AmpRelBreitWignerRes() 
{
}

void AmpRelBreitWignerRes::initialise() 
{
}

void AmpRelBreitWignerRes::setDecayMasses(double ma, double mb) {
  _ma = ma;
  _mb = mb;
}

double AmpRelBreitWignerRes::q0() const {
  double mapb = _ma + _mb;
  double mamb = _ma - _mb;
  
  return sqrt ( (_m0*_m0 - mapb*mapb) * (_m0*_m0 - mamb*mamb) ) / (2. * _m0 );
}

double AmpRelBreitWignerRes::q() const {
  double mapb = _ma + _mb;
  double mamb = _ma - _mb;
  
  return sqrt ( (_x*_x - mapb*mapb) * (_x*_x - mamb*mamb) ) / (2. * _x );
}


// compute part of the Blatt-Weisskopf barrier factor
//   BLprime = sqrt (F(q0)/F(q))
double AmpRelBreitWignerRes::F(double p) const {
  double retVal = 1;

  if (_spin == 0)
    retVal = 1;
  else if (_spin == 1) 
    retVal = 1 + p*p * _d*_d;
  else if (_spin == 2) {
    double z = p*p * _d*_d;
    retVal = (z-3.)*(z-3.) + 9*z;
  }
  return retVal;
}
    

// compute square of Blatt-Weisskopf barrier factor
double AmpRelBreitWignerRes::BLprime2() const {
  //  cout << q0() << " " << q() << "\t" << F(q0()) << " " << F(q()) << endl;
  return F(q0()) / F(q());
}


RooComplex AmpRelBreitWignerRes::evaluate() {

  // double m0 = p.Mass ();
  //    double Gamma0 = p.Width ();
  //    double q = p.q ();
  //    double q0 = p.q0 ();   
  //    double m = ~(p.get4P ());
  //    double GammaV;
  //    int l = p.Decay ()->L ();
 
  //    GammaV = Gamma0*(m0/m)*(q/q0)*(pow(F(l,q),2)/pow(F(l,q0),2));
  //    ret = (m0*Gamma0)/(m0*m0-m*m-complex<double>(0,1)*m0*GammaV);

  double Gamma0, GammaV;
  double m0 = Double_t(_m0);
  double m  = Double_t(_x);

  Gamma0 = double(_resWidth);
  GammaV = Gamma0 * (m0 / m) * pow(q() / q0(), 2.*_spin + 1.)  * BLprime2();
  
  RooComplex denom = RooComplex(m0*m0 - m*m, -m0 * GammaV);

  return (RooComplex(_m0 * Gamma0) / denom);
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
RooAbsArg* AmpRelBreitWignerRes::createFundamental(const char*) const {
  std::cout << "Unimplemented " << __PRETTY_FUNCTION__ << " called for " << this << std::endl
	    << "This might crash the program! " 
	    << std::endl;
  return 0;
}
