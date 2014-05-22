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

	//  std::vector<std::shared_ptr<AmpAbsDynamicalFunction> > _pdfList ;   //  List of component PDFs
	//  std::vector<std::shared_ptr<DoubleParameter> > _intList;    //  List of relative intensities
	//  std::vector<std::shared_ptr<DoubleParameter> > _phaseList;  //  List of relative phases
	//  std::vector<std::shared_ptr<AmpWigner> > _angList ;   //  List of component angular distributions
	//	std::cout<<"copy   "<<std::endl;

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

std::complex<double> AmpSumOfAmplitudes::getFirstBW(dataPoint& point) const
{
	//   RooComplex res;
	std::complex<double> res;
	//std::cout << "res = \t" << res.abs2() << std::endl;

	//std::cout << "PDFs: ";
	// for(unsigned int i=0; i<_pdfList.size(); i++){
	// double a = _intList[0]->GetValue();
	// double phi = _phaseList[0]->GetValue();
	// std::complex<double> eiphi(a*cos(phi),a*sin(phi));

	unsigned int twoJplusOne = (2*_pdfList[0]->getSpin()+1);
	//     res = res + (double)twoJplusOne * _pdfList[i]->evaluate() * eiphi;
	//twoJplusOne in included in evaluate(). We want to include this factor into the normalization of the amplitudes.
	res = (double)twoJplusOne * _pdfList[0]->evaluateAmp(point);
	// std::cout << _pdfList[i]->evaluate() << " ";
	//res = res + twoJplusOne * _pdfList[i]->evaluate() * eiphi * _angList[i]->evaluate();
	//}
	//std::cout << std::endl;

	return res;
}

std::complex<double> AmpSumOfAmplitudes::getFirstReso(dataPoint& point) const
{
	std::complex<double> res;
	//std::cout << "res = \t" << res.abs2() << std::endl;

	//std::cout << "PDFs: ";
	// for(unsigned int i=0; i<_pdfList.size(); i++){
	double a = _intList[0]->GetValue();
	double phi = _phaseList[0]->GetValue();
	std::complex<double> eiphi(a*cos(phi),a*sin(phi));

	double twoJplusOne = (2*_pdfList[0]->getSpin()+1);
	//     res = res + (double)twoJplusOne * _pdfList[i]->evaluate() * eiphi;
	res = twoJplusOne * _pdfList[0]->evaluate(point) * eiphi;
	// std::cout << _pdfList[i]->evaluate() << " ";
	//res = res + twoJplusOne * _pdfList[i]->evaluate() * eiphi * _angList[i]->evaluate();
	//}
	//std::cout << std::endl;

	return res;
}

std::complex<double> AmpSumOfAmplitudes::getFirstAmp(dataPoint& point) const
{
	std::complex<double> res;
	//std::cout << "res = \t" << res.abs2() << std::endl;

	//std::cout << "PDFs: ";
	for(unsigned int i=0; i<_pdfList.size(); i++){
		double a = _intList[i]->GetValue();
		double phi = _phaseList[i]->GetValue();
		std::complex<double> eiphi(a*cos(phi),a*sin(phi));

		double twoJplusOne = (2*_pdfList[i]->getSpin()+1);
		//     res = res + (double)twoJplusOne * _pdfList[i]->evaluate() * eiphi;
		res = res + twoJplusOne * _pdfList[i]->evaluate(point) * eiphi;
		// std::cout << _pdfList[i]->evaluate() << " ";
		//res = res + twoJplusOne * _pdfList[i]->evaluate() * eiphi * _angList[i]->evaluate();
	}
	//std::cout << std::endl;

	return res;
}

double AmpSumOfAmplitudes::evaluate(dataPoint& point) const
{
	std::complex<double> res;

	for(unsigned int i=0; i<_pdfList.size(); i++){
		double a = _intList[i]->GetValue();
		double phi = _phaseList[i]->GetValue();
		std::complex<double> eiphi(a*cos(phi),a*sin(phi));
		double twoJplusOne = (2*_pdfList[i]->getSpin()+1);
		double norm = _pdfList[i]->GetNormalization();

		//evaluate() = norm*evalAmp()*evalWignerD()
		res = res + twoJplusOne * _pdfList[i]->evaluate(point) * eiphi;
//		res = res + _pdfList[i]->evaluate(point) * eiphi;
	}

	return ( std::abs(res)*std::abs(res) );
}
double AmpSumOfAmplitudes::getUnormalizedFraction(std::string name){
	int id=-1;
	for(unsigned int i=0; i<_pdfList.size(); i++)
		if(_pdfList[i]->GetName()==name) id=i;
	if(id<0) return 0;
	return getUnormalizedFraction(id);
}
double AmpSumOfAmplitudes::getUnormalizedFraction(unsigned int id){
	double integral = _pdfList[id]->totalIntegral(); //integral over amplitude id: |amp*WignerD|^2
	double a = _intList[id]->GetValue();
	double phi = _phaseList[id]->GetValue();
	std::complex<double> eiphi(a*cos(phi),a*sin(phi));
	double coeff = std::abs(eiphi)*std::abs(eiphi);
	return ( coeff*integral );
}
double AmpSumOfAmplitudes::getTotalIntegral(std::string name){
	int id=-1;
	for(unsigned int i=0; i<_pdfList.size(); i++)
		if(_pdfList[i]->GetName()==name) id=i;
	if(id<0) return 0;
	return getTotalIntegral(id);
}
double AmpSumOfAmplitudes::getTotalIntegral(unsigned int id){
	std::cout<<"asdfasdfadfs"<<std::endl;
	unsigned int twoJplusOne = 2*_pdfList[id]->getSpin()+1;
	return ( _pdfList[id]->totalIntegral()*twoJplusOne );
//	return ( _pdfList[id]->totalIntegral() );
}

double AmpSumOfAmplitudes::evaluateSlice(dataPoint& point, std::complex<double>* reso, unsigned int nResos, unsigned int subSys=1) const
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
			res = res + reso[sys] * _angList[i]->evaluate(point);
		else
			res = res + _pdfList[i]->evaluate(point) * a * eiphi;
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

