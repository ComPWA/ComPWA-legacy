/*
 * momentumSystematics.cpp
 *
 *  Created on: Feb 27, 2015
 *      Author: weidenka
 */
#include <string>
#include <ostream>

#include "Core/Logging.hpp"
#include "Core/TableFormater.hpp"
#include "DataReader/CorrectionTable.hpp"

void CorrectionTable::Print() const{
	if(sys.size() == 0 && antiSys.size() == 0) return; //don't print if empty
	std::stringstream out;
	TableFormater table(&out);
	table.addColumn("Bin",12);//add empty first column
	table.addColumn("Particle",20);//global correlation coefficient
	table.addColumn("anti-Particle",20);//global correlation coefficient
	table.addColumn("Neutral/Average",20);//global correlation coefficient
	table.header();
	for(unsigned int i=0; i<bins.size(); i++){
		std::stringstream strBin; strBin << bins.at(i).first << " - " << bins.at(i).second;
		table << strBin.str();
		std::stringstream strCor;
		strCor << GetValue(+1,bins.at(i).second);
//				<< " +- " <<GetError(+1,bins.at(i).second);
		table << strCor.str();
		std::stringstream strAntiCor;
		strAntiCor << GetValue(-1,bins.at(i).second);
//				<< " +- " <<GetError(-1,bins.at(i).second);
		table << strAntiCor.str();
		std::stringstream strAvg;
		strAvg << GetValue(0,bins.at(i).second);
//				<< " +- " <<GetError(0,bins.at(i).second);
		table << strAvg.str();
	}
	table.footer();
	BOOST_LOG_TRIVIAL(info) << "CorrectionTable::Print() | "<<title << std::endl << out.str();

}
void CorrectionTable::SetSystematics(std::vector<double> b, std::vector<double> bError ){
	sys = b;
	if(bError.size() == b.size())//is error vector valid?
		sysError = bError;
	else //otherwise set errors to 0
		sysError = std::vector<double>(bins.size(),0.);
	if(!antiSys.size()) {//if anti-particle systematics is empty set is to the same values
		antiSys = b;
		antiSysError = bError;
	}
	return;
}
void CorrectionTable::SetAntiSystematics(std::vector<double> b, std::vector<double> bError ){
	antiSys = b;
	if(bError.size() == b.size())//is error vector valid?
		antiSysError = bError;
	else //otherwise set errors to 0
		antiSysError = std::vector<double>(bins.size(),0.);
	return;
}
void CorrectionTable::SetSystematicsError(std::vector<double> bError){
	sysError = bError;
	return;
}
void CorrectionTable::SetAntiSystematicsError(std::vector<double> bError){
	antiSysError = bError;
	return;
}
void CorrectionTable::AddInv(double binMin, double binMax, double s, double sError,
		double antiS, double antiSerror){
	Add(binMin, binMax, inverse(s), inverseError(s, sError),
			inverse(antiS), inverseError(antiS,antiSerror) );
}
void CorrectionTable::Add(double binMin, double binMax, double s, double sError,
		double antiS, double antiSerror){
	addBin(binMin, binMax);
	sys.push_back(s);
	sysError.push_back(sError);
	if(antiS==-999){//no efficiency is given for anti-particles
		antiSys.push_back(s);
		antiSysError.push_back(sError);
	} else {
		antiSys.push_back(antiS);
		antiSysError.push_back(antiSerror);
	}
}
double CorrectionTable::GetInvValue(int charge, double momentum) const{
	return inverse(GetValue(charge, momentum));
}
double CorrectionTable::GetInvError(int charge, double momentum) const{
	return inverseError(GetValue(charge, momentum),GetError(charge,momentum));
}
double CorrectionTable::GetValue(int charge, double momentum) const{
	int binNumber = findBin(momentum);
	if(binNumber<0)
		throw std::runtime_error("momentumSys::GetValue() no bin found for momentum value "
				+ std::to_string((long double)momentum)+"!");
	if(charge>0)//D0->K0bar K+K-
		return sys[binNumber];
	else if(charge<0)//D0bar->K0 K+K-
		return antiSys[binNumber];
	else if(charge==0){//charge unknown -> average value
		return (antiSys[binNumber]+sys[binNumber])/2;
	}
	return 0.;
}
double CorrectionTable::GetError(int charge, double momentum) const{
	int binNumber = findBin(momentum);
	if(binNumber<0)
		throw std::runtime_error("momentumSys::GetError() no bin found for momentum value "
				+ std::to_string((long double)momentum)+"!");
	if(charge>0)//D0->K0bar K+K-
		return sysError[binNumber];
	else if(charge<0)//D0->K0bar K+K-
		return antiSysError[binNumber];
	else if(charge==0){//charge unknown -> average value
		return (0.5*std::sqrt(
				antiSysError[binNumber]*antiSysError[binNumber]
													 +sysError[binNumber]*sysError[binNumber])
		);
	}
	return 0.;
}

void CorrectionTable::AddToTotalError(int charge, double momentum){
	try{
		totalSys.push_back( GetValue(charge, momentum) );
		totalSysError.push_back( GetError(charge, momentum) );
	} catch ( const std::exception &e) {
		BOOST_LOG_TRIVIAL(info) << "momentumSys::AddToTotalError()  can't add "
				"systematics for momentum "<<momentum<<" to total systematics. We skip that track!";
	}
}
double CorrectionTable::GetTotalSystematics(bool useArithmeticMean) {
	double mean=0;
	double sumSigma = 0;
	//use arithmetic mean
	if(useArithmeticMean){
		for(unsigned int i=0; i<totalSys.size(); i++)
			mean+=totalSys.at(i);
		return mean/totalSys.size();
	} else {
		//use weighted mean
		for(unsigned int i=0; i<totalSys.size(); i++){
			double tmp = 1/( totalSysError.at(i)*totalSysError.at(i) );
			mean+=tmp*totalSys.at(i);
			sumSigma+=tmp;
		}
		mean/=sumSigma;
		return mean;
	}
}
double CorrectionTable::GetTotalSystematicsError(bool useArithmeticMean) {
	double sumSigma = 0; //inverse sum over uncertainties
	double mean = GetTotalSystematics(useArithmeticMean);
	if(useArithmeticMean){
		for(unsigned int i=0; i<totalSys.size(); i++)
			sumSigma += (totalSys.at(i)-mean)*(totalSys.at(i)-mean);
		return sumSigma = sqrt(sumSigma/totalSys.size());
	} else {
		for(unsigned int i=0; i<totalSysError.size(); i++)
			sumSigma+= 1/( totalSysError.at(i)*totalSysError.at(i) );
		return sqrt(1/sumSigma);
	}
}

//! invert e_mc/e_data-1 to e_data/e_mc-1
double CorrectionTable::inverse(double x){
	if(x==-999 || x==-1) return -999;
	return 1/(x+1)-1;
}
//! Calculate error for inversion e_mc/e_data-1 to e_data/e_mc-1
double CorrectionTable::inverseError(double x, double xErr){
	if(x==-999 || x==-1) return -999;
	return std::abs(-1/((x+1)*(x+1)) * xErr);
}
int CorrectionTable::findBin(double momentum) const {
	for(unsigned int i=0; i<bins.size(); i++){
		if( momentum <= bins.at(i).second && momentum > bins.at(i).first){
			return i;
		}
	}
	return -1;//no bin found
}
//! check if there is a bin defined that overlaps with [min,max]
int CorrectionTable::findBin(double min, double max) const {
	for(unsigned int i=0; i<bins.size(); i++){
		if( (min < bins.at(i).second && min >= bins.at(i).first) ||
				( max <= bins.at(i).second	&& max > bins.at(i).first ) ) return i;
	}
	return -1;//bin doesn't exist
}
void CorrectionTable::addBin(double min, double max){
	if( findBin(min,max) < 0 ) //is a bin already defined?
		bins.push_back(std::make_pair(min, max));
	else throw std::runtime_error("CorrectionTable::addBin() in range ["
			+std::to_string((long double)min)+","+std::to_string((long double)max)+"] a bin is already defined!");
}
bool CorrectionTable::check() const{
	int s = bins.size();
	if(sys.size() != s || sysError.size() != s || antiSys.size() != s || antiSysError.size() != s)
		return 0;
	return 1;
}
