//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
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


AmpFlatteRes::AmpFlatteRes(const AmpFlatteRes& other, const char* newname) :
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

AmpFlatteRes::AmpFlatteRes(const AmpFlatteRes& other) :
  AmpAbsDynamicalFunction(other.GetName(), other.GetTitle()),
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

AmpFlatteRes::~AmpFlatteRes() 
{
}

void AmpFlatteRes::initialise() 
{
}

void AmpFlatteRes::setDecayMasses(double ma, double mb) {
  _ma = ma;
  _mb = mb;
}

void AmpFlatteRes::setBarrierMass(double mBarA, double mBarB) {
  _mBarA = mBarA;
  _mBarB = mBarB;
}

double AmpFlatteRes::q0() const {
  double mapb = _ma + _mb;
  double mamb = _ma - _mb;
  
  return sqrt ( (_m0*_m0 - mapb*mapb) * (_m0*_m0 - mamb*mamb) ) / (2. * _m0 );
}

double AmpFlatteRes::q() const {
  double mapb = _ma + _mb;
  double mamb = _ma - _mb;
  
  return sqrt ( (_x*_x - mapb*mapb) * (_x*_x - mamb*mamb) ) / (2. * _x );
}

double AmpFlatteRes::q0(double ma, double mb) const {
  double mapb = ma + mb;
  double mamb = ma - mb;
  
  return sqrt ( (fabs(_m0*_m0 - mapb*mapb)) * (_m0*_m0 - mamb*mamb) ) / (2. * _m0 );
}

double AmpFlatteRes::q(double ma, double mb) const {
  double mapb = ma + mb;
  double mamb = ma - mb;
  
  return sqrt ( (_x*_x - mapb*mapb) * (_x*_x - mamb*mamb) ) / (2. * _x );
}


// compute part of the Blatt-Weisskopf barrier factor
//   BLprime = sqrt (F(q0)/F(q))
double AmpFlatteRes::F(double p) const {
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
double AmpFlatteRes::BLprime2() const {
  //  cout << q0() << " " << q() << "\t" << F(q0()) << " " << F(q()) << endl;
  return F(q0()) / F(q());
}

/*double AmpFlatteRes::phspTerm(double ma, double mb) const {
  double mapb = ma + mb;
  double mamb = ma - mb;
  
  return sqrt ( (1.-(mamb*mamb/Double_t(_x)/Double_t(_x)))*(1.-(mapb*mapb/Double_t(_x)/Double_t(_x))) );
}*/

RooComplex AmpFlatteRes::evaluate() const {

  double m0 = Double_t(_m0);
  double m  = Double_t(_x);

  double Gamma0 = double(_resWidth);
  //double GammaV = Gamma0 * (m0 / m) * pow(q() / q0(), 2.*_spin + 1.)  * BLprime2();
  double gam1=0.6325, gam2=0.7746; //gam1=1./sqrt(2./5.), gam2=1./sqrt(3./5.); //gam1^2+gam2^2=1; TODO: set channel dependent
  double p1 = q0(_ma,_mb);
  double p2 = q0(_mBarA,_mBarB);
  double g1 = gam1*sqrt(m0*Gamma0);
  double g2 = gam2*sqrt(m0*Gamma0);
  
  RooComplex denom = RooComplex(m0*m0 - m*m, -p1*g1*g1-p2*g2*g2);

  RooComplex result = (RooComplex(_m0 * Gamma0) / denom); //sqrt(Gamma0 * GammaI) ?
  if(result.re()!=result.re()) {std::cout << "OH OH RE" << std::endl; return 0;}
  if(result.im()!=result.im()) {std::cout << "OH OH IM" << std::endl; return 0;}
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
