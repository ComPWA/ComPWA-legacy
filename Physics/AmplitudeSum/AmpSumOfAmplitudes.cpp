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
#include "Core/Parameter.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

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

void AmpSumOfAmplitudes::addBW(std::shared_ptr<AmpAbsDynamicalFunction> theRes , std::shared_ptr<DoubleParameter> r, std::shared_ptr<DoubleParameter> phi, std::shared_ptr<AmpWigner2> theAng) {
	_pdfList.push_back(theRes);
	_intList.push_back(r);
	_phaseList.push_back(phi);
}

void AmpSumOfAmplitudes::addBW(std::shared_ptr<AmpAbsDynamicalFunction> theRes , std::shared_ptr<DoubleParameter> r, std::shared_ptr<DoubleParameter> phi) {
	_pdfList.push_back(theRes);
	_intList.push_back(r);
	_phaseList.push_back(phi);
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

double AmpSumOfAmplitudes::evaluate(const dataPoint& point) const
{
	std::complex<double> res;
	for(unsigned int i=0; i<_pdfList.size(); i++){
		double a = _intList[i]->GetValue();
		double phi = _phaseList[i]->GetValue();
		std::complex<double> eiphi(a*cos(phi),a*sin(phi));
		res = res + _pdfList[i]->evaluate(point) * eiphi;//adding factor 2J+1 to AmpWigner2 -> consistency with tree
	}
	return ( std::norm(res) ); //return |A|^2
}

double AmpSumOfAmplitudes::getSpin(unsigned int id){
	return _pdfList[id]->getSpin();
}

double AmpSumOfAmplitudes::getAmpStrength(unsigned int id){
	double integral = getAmpIntegral(id);
	double mag = _intList[id]->GetValue();
	return ( mag*mag*integral );
}

double AmpSumOfAmplitudes::getAmpIntegral(unsigned int id){
//	return ( _pdfList[id]->totalIntegral() ); //2J+1 is already in AmpWigner2
	return ( 2*_pdfList[id]->getSpin()+1); //assume that amplitude is normalized, save some cpu time
}

int AmpSumOfAmplitudes::getAmpId(std::string name) {
	for(unsigned int i=0; i< _pdfList.size(); i++)
		if(_pdfList.at(i)->GetName()==name) return i;
	BOOST_LOG_TRIVIAL(error) << "AmpSumOfAmplitudes::getAmpId() | amplitude with name "<<name<<" not found in list!";
	return -999;
}

double AmpSumOfAmplitudes::evaluateSlice(dataPoint& point, std::complex<double>* reso, unsigned int nResos, unsigned int subSys=5, double N=1., unsigned int nF0=3, unsigned int nF2=2) const
{
	// ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
	std::complex<double> res;

	int sys=0;

	//std::cout<< "DEBUG Point " << point.getVal("m23sq") << " " << point.getVal("m13sq") << " " << reso[0] << " " << reso[1] << std::endl;

	bool used[nResos];
	for(unsigned int i=0; i<nResos; i++) used[i]=false;

	for(unsigned int i=0; i<_pdfList.size(); i++){
		double a = _intList[i]->GetValue();
		double phi = _phaseList[i]->GetValue();
		std::complex<double> eiphi (a * cos(phi), a * sin(phi));
		double twoJplusOne = (2*_pdfList[i]->getSpin()+1);
		if(i<nF0) sys = 0; //TODO: way better!!!
		else if(i<(nF2+nF0)) sys = 1;
		else sys = 999;
		//sys = itReso;

		//std::cout<< "DEBUG Reso " << i << " spinsys " << sys << " " << _pdfList[i]->isSubSys(subSys) << " use " << used[sys] << std::endl;

		if(dynamic_cast<AmpKinematics*>(&*(_pdfList[i]))->isSubSys(subSys)){
		  if(!used[sys]){
			res = res + reso[sys] * dynamic_cast<AmpKinematics*>(&*(_pdfList[i]))->evaluateWignerD(point);
		//if(_pdfList[i]->isSubSys(subSys)){
		//  if(!used[sys] && sys<999){ //add terms only ones
		//	res = res + N * twoJplusOne * reso[sys] * _angList[i]->evaluate(point);
			used[sys]=true;
		  }
		}else{
			res = res + N * twoJplusOne * _pdfList[i]->evaluate(point) * eiphi;
			//std::cout<< " Omega " << i << std::endl;
		}
		//res = res + twoJplusOne * _pdfList[i]->evaluate(point) * eiphi;
		//itReso++;
	}

	return (std::abs(res)*std::abs(res) );
}

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */
