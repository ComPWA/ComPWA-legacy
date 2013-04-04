//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpWigner.hpp"
#include "RooAbsReal.h"
#include "RooRealVar.h"

#include "qft++.h"

AmpWigner::AmpWigner():
  RooAbsArg("none","")
{
  toEvaluate=false;
}


AmpWigner::AmpWigner(const char *name, const char *title,
		       RooAbsReal& m12, ///  
                       RooAbsReal& m23, ///  
                       RooAbsReal& m13, ///  
                       const int subSysFlag, ///  
		       RooAbsReal& inSpin, RooAbsReal& outSpin1,
		       RooAbsReal& outSpin2) :
  RooAbsArg(name,title),
  _m12("m12", "Observable", this, m12),
  _m23("m23", "Observable", this, m23),
  _m13("m13", "Observable", this, m13),
  _subSysFlag(subSysFlag),
  _inSpin("inSpin", "inSpin", this, inSpin),
  _outSpin1("outSpin1", "outSpin", this, outSpin1),
  _outSpin2("outSpin2", "outSpin", this, outSpin2)
{
  toEvaluate=true;
  initialise();

  _M = 3.096916;  //TODO: Setter!
  _m1 = 0.;
  _m2 = 0.139570;
  _m3 = 0.139570;
}


AmpWigner::AmpWigner(const AmpWigner& other, const char* newname) :
  RooAbsArg(other, newname),
  _m12("m12", "Observable", this, other._m12),
  _m23("m23", "Observable", this, other._m23),
  _m13("m13", "Observable", this, other._m13),
  _subSysFlag(other._subSysFlag),
  _inSpin("inSpin", this, other._inSpin),
  _outSpin1("outSpin1", this, other._outSpin1),
  _outSpin2("outSpin2", this, other._outSpin2),
  toEvaluate(other.toEvaluate)
{
  initialise();

  _M = 3.096916;  //TODO: Setter!
  _m1 = 0.;
  _m2 = 0.139570;
  _m3 = 0.139570;
}

AmpWigner::~AmpWigner() 
{
}

void AmpWigner::initialise() 
{
}    

double AmpWigner::lambda(double x, double y, double z)const{
  return x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z;
}

Double_t AmpWigner::s2min(Double_t s1, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const
{
	Double_t s      = m0*m0;
	Double_t lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );
	 
	Double_t result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) - lamterm )/(2.*s1);
	
	return result;
}

Double_t AmpWigner::s2max(Double_t s1, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const
{
	Double_t s      = m0*m0;
	Double_t lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m2*m2, m3*m3) );
	 
	Double_t result  = m1*m1 + m3*m3 + ( (s-s1-m1*m1)*(s1-m2*m2+m3*m3) + lamterm )/(2.*s1);
	
	return result;
}

Double_t AmpWigner::s3min(Double_t s1, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const
{
	Double_t s      = m0*m0;
	Double_t lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );
	 
	Double_t result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m2*m2) - lamterm )/(2.*s1);
	
	return result;
}

Double_t AmpWigner::s3max(Double_t s1, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const
{
	Double_t s      = m0*m0;
	Double_t lamterm = sqrt( lambda(s1,s,m1*m1) ) * sqrt( lambda(s1, m3*m3, m1*m1) );
	 
	Double_t result  = m1*m1 + m2*m2 + ( (s-s1-m1*m1)*(s1-m1*m1+m3*m3) + lamterm )/(2.*s1);
	
	return result;
}

Double_t AmpWigner::s1min(Double_t s2, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const
{
	Double_t s      = m0*m0;
	Double_t lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m3*m3, m1*m1) );
	 
	Double_t result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) - lamterm )/(2.*s2);
	
	return result;
}

Double_t AmpWigner::s1max(Double_t s2, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const
{
	Double_t s      = m0*m0;
	Double_t lamterm = sqrt( lambda(s2,s,m2*m2) ) * sqrt( lambda(s2, m1*m1, m3*m3) );
	 
	Double_t result  = m2*m2 + m3*m3 + ( (s-s2-m2*m2)*(s2-m1*m1+m3*m3) + lamterm )/(2.*s2);
	
	return result;
}

double AmpWigner::evaluate() const {
  if(!toEvaluate) return 1.;

  double locmin_sq, locmax_sq, beta;

  switch(_subSysFlag){
    case 1:{ //reso in m23
      locmin_sq = s2min(_m23*_m23,_M,_m1,_m2,_m3); 
      locmax_sq = s2max(_m23*_m23,_M,_m1,_m2,_m3);
      beta=acos((2.*_m13*_m13-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
      break;
    }
    case 2:{ //reso in m13
      locmin_sq = s1min(_m13*_m13,_M,_m1,_m2,_m3); 
      locmax_sq = s1max(_m13*_m13,_M,_m1,_m2,_m3);
      beta=acos((2.*_m23*_m23-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
      break;
    }
    case 3:{ //reso in m12
      //return 1;
      locmin_sq = s1min(_m12*_m12,_M,_m1,_m3,_m2); 
      locmax_sq = s1max(_m12*_m12,_M,_m1,_m3,_m2);
      beta=acos((2.*_m23*_m23-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
      if(beta!=beta) return 1.;
      break;
    }
  }

  //double locmin_sq = s2min(_y*_y), locmax_sq = s2max(_y*_y);
  //if( _x*_x>locmax_sq || _x*_x<locmin_sq )
  //  return 0.;

  Spin j(_inSpin), m(_outSpin1), n(_outSpin2);
  return Wigner_d(j,m,n,beta);
}
                       
////  implement stubs for virtual functions
////  to satisfy interface of RooAbsArg.
////  not sure which ones are really needed.

TObject*  AmpWigner::clone (const char *newname) const 
{
  //  std::cout << "Unimplemented AmpRelBreitWignerRes::clone called for " << this << endl
  //	    << "This will crash the program! " 
  //	    << endl;

  //  cout << __PRETTY_FUNCTION__ << "(" << ((newname)?newname:"(null)") << ")" << endl;
  return new AmpWigner(*this, newname);
}

Bool_t AmpWigner::readFromStream(std::istream&, Bool_t, Bool_t) {
  cout << __PRETTY_FUNCTION__ << endl;
  return false;
}

void AmpWigner::writeToStream(std::ostream&, Bool_t) const {
  cout << __PRETTY_FUNCTION__ << endl;
 ; 
}

Bool_t AmpWigner::operator==(const RooAbsArg&) {
  cout << __PRETTY_FUNCTION__ << endl;
  return false;
}

void AmpWigner::syncCache(const RooArgSet*)  {
  cout << __PRETTY_FUNCTION__ << endl;
  ;
}

void AmpWigner::copyCache(const RooAbsArg*, Bool_t, Bool_t) {
  cout << __PRETTY_FUNCTION__ << endl;
  ; 
}

void AmpWigner::attachToTree(TTree&, Int_t) {
  cout << __PRETTY_FUNCTION__ << endl;
  ;
}

void AmpWigner::attachToVStore(RooVectorDataStore&) {
  cout << __PRETTY_FUNCTION__ << endl;
  ;
}

void AmpWigner::setTreeBranchStatus(TTree&, Bool_t) {
  cout << __PRETTY_FUNCTION__ << endl;
  ;
}

void AmpWigner::fillTreeBranch(TTree&) {
  cout << __PRETTY_FUNCTION__ << endl;
  ;
}
RooAbsArg* AmpWigner::createFundamental(const char*) const {
  std::cout << "Unimplemented " << __PRETTY_FUNCTION__ << " called for " << this << endl
	    << "This might crash the program! " 
	    << endl;
  return 0;
}
