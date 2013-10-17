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
		UInt_t spin, UInt_t m,  UInt_t n,UInt_t subSys) :
						  RooAbsArg(name,title),
						  //				  _m13("m13", "Observable", 0.0),
						  //				  _m23("m23", "Observable", 0.0),
						  _inSpin(spin),
						  _outSpin1(m),
						  _outSpin2(n),
						  _subSys(subSys)
{
	toEvaluate=true;
	initialise();
}


AmpWigner::AmpWigner(const char *name, const char *title,
		RooAbsReal& m13, ///
		RooAbsReal& m23, ///
		UInt_t spin, UInt_t m,  UInt_t n,UInt_t subSys) :
						  RooAbsArg(name,title),
						  _m13("m13", "Observable", this, m13),
						  _m23("m23", "Observable", this, m23),
						  _inSpin(spin),
						  _outSpin1(m),
						  _outSpin2(n),
						  _subSys(subSys)
{
	toEvaluate=true;
	initialise();
}


AmpWigner::AmpWigner(const AmpWigner& other, const char* newname) :
						  RooAbsArg(other, newname),
						  _m13("m13", "Observable", this, other._m13),
						  _m23("m23", "Observable", this, other._m23),
						  _inSpin(other._inSpin),
						  _outSpin1(other._outSpin1),
						  _outSpin2(other._outSpin2),
						  _subSys(other._subSys),
						  toEvaluate(other.toEvaluate)
{
	initialise();
}

AmpWigner::~AmpWigner() 
{
}

void AmpWigner::initialise() 
{
	static dataPoint* point = dataPoint::instance();
//	cout<<"================== "<<point->DPKin.m3<<endl;
	_M=point->DPKin.M;
	if(_subSys==5){
		_m1=point->DPKin.m1;
		_m2=point->DPKin.m2;
		_m3=point->DPKin.m3;}
	if(_subSys==4){
		_m1=point->DPKin.m2;
		_m2=point->DPKin.m3;
		_m3=point->DPKin.m1;}
	if(_subSys==3){
		_m1=point->DPKin.m3;
		_m2=point->DPKin.m2;
		_m3=point->DPKin.m1;}
}    

void AmpWigner::setDecayMasses(double m1, double m2, double m3, double M){
	_M=M; _m1=m3; _m2=m2; _m3=m1;
	return;
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

//	cout<<"==== "<<_subSys<< " " <<_m1<< " "<<_m2<<" " <<_m3<<endl;
	double invM1 = dataPoint::instance()->getM(_subSys);
	/*
	 * For WignerD functions we need one more invariant mass:
	 */
	int mod=0;
	if(_subSys==5) mod=4; //5->3 work also without beta=nan, what is correct?
	if(_subSys==4) mod=3;
	if(_subSys==3) mod=4;
	double invM2 = dataPoint::instance()->getM(mod);

	dataPoint* point = dataPoint::instance();
	//	cout<<point->getM(3)<<" " <<point->getM(4)<< " " << point->getM(5)<<endl;
	//	double locmin_sq2 = s2min(_m23*_m23,_M,_m1,_m2,_m3);
	//	double locmax_sq2 = s2max(_m23*_m23,_M,_m1,_m2,_m3);
	//	double beta2=acos((2.*_m13*_m13-locmax_sq2-locmin_sq2)/(locmax_sq2-locmin_sq2));

	locmin_sq = s2min(invM1*invM1,_M,_m1,_m2,_m3);
	locmax_sq = s2max(invM1*invM1,_M,_m1,_m2,_m3);
	beta=acos((2.*invM2*invM2-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
	//	cout<<"==== "<<_m23<< " "<<_m13<< " "<<invM1<<" " <<invM2<<endl;
	//	cout<< "wwww " <<_subSys<< " "<<mod<<" " <<beta<< " "<<beta2<<endl;
	//	if(beta!=beta){
	//		std::cout<<beta<<std::endl;
	//		std::cout<<_M<< " "<<_m1<<" " <<_m2<<" " <<_m3<<std::endl;
	//		std::cout<<(2.*_m13*_m13-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq)<<std::endl;
	//		std::cout<<_m13<< " " <<_m23<<" " <<locmin_sq << " "<<locmax_sq<<std::endl;
	//	}
	//double locmin_sq = s2min(_y*_y), locmax_sq = s2max(_y*_y);
	//if( _x*_x>locmax_sq || _x*_x<locmin_sq )
	//  return 0.;

	Spin j((double) _inSpin), m((double) _outSpin1), n((double)_outSpin2);
	//  cout<<Wigner_d(j,m,n,beta)<<endl;
	double result = Wigner_d(j,m,n,beta);
	if( ( result!=result ) || (beta!=beta)) {
		std::cout<< "NAN! J="<< _inSpin<<" M="<<_outSpin1<<" N="<<_outSpin2<<" beta="<<beta<<std::endl;
		std::cout<< "subSys: "<<_subSys<<" ("<<mod<<") "<<invM1 << " " <<invM2<< " cos(beta)="<<(2.*invM2*invM2-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq)<<std::endl;
		return 0;
	}
	//	cout<<"result: "<<result<<endl;
	return result;
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

Bool_t AmpWigner::isIdentical(const RooAbsArg&, Bool_t){
	cout << __PRETTY_FUNCTION__ << endl;
	;
}
