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
#include <cmath> 

//#include "Riostream.h"
//#include "TMath.h"

//#include "RooAbsReal.h"
#include "Core/Parameter.hpp"
//#include "RooAbsCategory.h"
//#include "RooLinkedListIter.h"
//#include "boost/function.hpp"
//#include "boost/bind.hpp"
#include <functional>

#include "qft++.h"

#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"

//ClassImp(AmpSumOfAmplitudes);

AmpSumOfAmplitudes::AmpSumOfAmplitudes()
{

}

 AmpSumOfAmplitudes::AmpSumOfAmplitudes(const char *name)
 { 

 } 


 AmpSumOfAmplitudes::AmpSumOfAmplitudes(const AmpSumOfAmplitudes& other, const char* name)
 { 

 } 

 AmpSumOfAmplitudes::~AmpSumOfAmplitudes(){
   //something TODO?
 }

void AmpSumOfAmplitudes::addBW(std::shared_ptr<AmpAbsDynamicalFunction> theRes , std::shared_ptr<DoubleParameter> r, std::shared_ptr<DoubleParameter> phi, std::shared_ptr<AmpWigner> theAng) {
  _pdfList.push_back(theRes);
  _intList.push_back(r);
  _phaseList.push_back(phi);
  _angList.push_back(theAng);
}

void AmpSumOfAmplitudes::addBW(std::shared_ptr<AmpAbsDynamicalFunction> theRes , std::shared_ptr<DoubleParameter> r, std::shared_ptr<DoubleParameter> phi) {
  _pdfList.push_back(theRes);
  _intList.push_back(r);
  _phaseList.push_back(phi);
  _angList.push_back(std::shared_ptr<AmpWigner>(new AmpWigner()));
}

double AmpSumOfAmplitudes::integrate(AmpRelBreitWignerRes* amp){

	using namespace std::placeholders;

	double x[2]; x[0]=1.3; x[1]=1.3;
//	auto d = std::bind(&AmpRelBreitWignerRes::eval,*amp,_1,_2,_3);
//	cout<<d(x,2,0)<<std::endl;
//	auto d = std::bind(&AmpRelBreitWignerRes::eval,*amp,x,2,0);
//	cout<<d()<<std::endl;

//	 boost::function<int (int)> f;
//	X x;
//	f = std::bind1st(
//	      std::mem_fun(&X::foo), &x);
//	f(5); // Call x.foo(5)

//	std::cout<<" " <<sp(this,x,2,0)<<std::endl;
//	std::cout<<"====12=== "<<bb(x,2,0)<<std::endl;
//		gsl_monte_function G = {sp,2,0};
}
 double AmpSumOfAmplitudes::evaluate() const
 { 
//   RooComplex res;
   complex<double> res;
   //std::cout << "res = \t" << res.abs2() << std::endl;

   //std::cout << "PDFs: ";
   for(unsigned int i=0; i<_pdfList.size(); i++){
     double a = _intList[i]->GetValue();
     double phi = _phaseList[i]->GetValue();
     std::complex<double> eiphi(a*cos(phi),a*sin(phi));

     std::complex<double> twoJplusOne(2*_pdfList[i]->getSpin()+1);
     res = res + twoJplusOne * _pdfList[i]->evaluate() * eiphi;
    // std::cout << _pdfList[i]->evaluate() << " ";
//res = res + twoJplusOne * _pdfList[i]->evaluate() * eiphi * _angList[i]->evaluate();
   }
   //std::cout << std::endl;

   return ( std::abs(res)*std::abs(res) );
 } 

 double AmpSumOfAmplitudes::evaluateSlice(std::complex<double>* reso, unsigned int nResos, unsigned int subSys=1) const
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   std::complex<double> res;

   int itReso=0, sys=0;

   for(unsigned int i=0; i<_pdfList.size(); i++){
     double a = _intList[i]->GetValue();
     double phi = _phaseList[i]->GetValue();
     std::complex<double> eiphi (cos(phi), sin(phi));
     if(itReso<3) sys = 0; //TODO: way better!!!
     else if(itReso<5) sys = 1;
     //else sys = 2;
     //sys = itReso;

     if(_pdfList[i]->isSubSys(subSys))
       res = res + reso[sys] * _angList[i]->evaluate();
     else
       res = res + _pdfList[i]->evaluate() * a * eiphi;
 //res = res + _pdfList[i]->evaluate() * a * eiphi * _angList[i]->evaluate();
     itReso++;
   }

   return (std::abs(res)*std::abs(res) );
 } 

 /*double AmpSumOfAmplitudes::evaluatePhi() const
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   RooComplex res;

   AmpRelBreitWignerRes *pdf;
   DoubleParameter *theInt;
   DoubleParameter *thePhase;
   AmpWigner *ang;

   _pdfIter->Reset();
   _intIter->Reset();
   _phaseIter->Reset();
   _angIter->Reset();

   //   TIterator* _pdfIter = _pdfList.createIterator() ;
   //   AmpRelBreitWignerRes *pdf;


   while((pdf      = (AmpRelBreitWignerRes*)_pdfIter->Next()) &&
	 (theInt   = (DoubleParameter*)_intIter->Next())        &&
	 (thePhase = (DoubleParameter*)_phaseIter->Next())      &&
         (ang      = (AmpWigner*)_angIter->Next())  ) {
     double a = theInt->GetValue();
     double phi = thePhase->GetValue();
     RooComplex eiphi (cos(phi), sin(phi));

     res = res + pdf->evaluate() * a * eiphi * ang->evaluate();
   }

   return atan2(res.im(),res.re()); 
 } */

