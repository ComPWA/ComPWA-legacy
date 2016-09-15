/*
 * FitResult.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: weidenka
 */
#include "Core/FitResult.hpp"

namespace ComPWA {

void FitResult::writeText(std::string filename){
	std::ofstream myfile;
	myfile.open(filename, std::ios::app);
	genOutput(myfile);
	myfile.close();
	return;
};

void FitResult::writeSimpleText(std::string filename){
	std::ofstream myfile;
	myfile.open(filename);
	genSimpleOutput(myfile);
	myfile.close();
	return;
};

double FitResult::shiftAngle(double v){
	double originalVal = v;
	double val = originalVal;
	double pi = ComPWA::PhysConst::Instance().findConstant("Pi").value_;
	while(val> pi) val-=2*pi;
	while(val< -pi ) val+=2*pi;
	if(val!=originalVal)
		BOOST_LOG_TRIVIAL(info) << "shiftAngle(): shifting parameter from "<<originalVal<< " to "<<val<<"!";
	return val;
}

void FitResult::genSimpleOutput(std::ostream& out)
{
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> outPar = finalParameters.GetDoubleParameter(o);
		out<<outPar->GetValue()<<" "<<outPar->GetError()<<" ";
	}
	out<<"\n";

	return;
}

void FitResult::setFinalParameters(ParameterList finPars)
{
	finalParameters.DeepCopy(finPars);
}

void FitResult::setUseCorrelatedErrors(int nSets)
{
	if(nSets < 0) nSetsFractionError = 0;
	nSetsFractionError = nSets;
	return;
}

void FitResult::print(std::string opt)
{
	std::stringstream s;
	genOutput(s,opt);
	std::string str = s.str();
	BOOST_LOG_TRIVIAL(info) << str;
}

void FitResult::printFitParameters(TableFormater* tableResult)
{
	bool printTrue=0, printInitial=0;
	if(trueParameters.GetNParameter()) printTrue=1;
	if(initialParameters.GetNParameter()) printInitial=1;

	//Column width for parameter with symmetric error
	unsigned int parErrorWidth = 22;

	//Do we have a parameter with assymetric errors?
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++)
		if(finalParameters.GetDoubleParameter(o)->GetErrorType()==ErrorType::ASYM)
			parErrorWidth=33;

	tableResult->addColumn("Nr");
	tableResult->addColumn("Name",15);
	if(printInitial) tableResult->addColumn("Initial Value",parErrorWidth);
	tableResult->addColumn("Final Value",parErrorWidth);
	if(printTrue) tableResult->addColumn("True Value",10);
	if(printTrue) tableResult->addColumn("Deviation",9);
	tableResult->header();

	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> iniPar, outPar, truePar;
		try{
			outPar = finalParameters.GetDoubleParameter(o);
		} catch (BadParameter& bad){
			BOOST_LOG_TRIVIAL(error) << "FitResult::printFitParameters() | "
					"can't access parameter of final parameter list!";
			throw;
		}
		if(printInitial){
			try{
				iniPar = initialParameters.GetDoubleParameter(outPar->GetName());
			} catch (BadParameter& bad){
				iniPar = 0;
			}
		}
		if(printTrue){
			try{
				truePar = trueParameters.GetDoubleParameter(outPar->GetName());
			} catch (BadParameter& bad){
				iniPar = 0;
			}
		}

		ErrorType errorType = outPar->GetErrorType();
		bool isFixed = outPar->IsFixed();

		// Is parameter an angle?
		bool isAngle=0;
		if(outPar->GetName().find("phase")!=string::npos) isAngle=1;
		// ... then shift the value to the domain [-pi;pi]
		if(isAngle && !isFixed) {
			outPar->SetValue( shiftAngle(outPar->GetValue()) );
			if(printInitial) iniPar->SetValue( shiftAngle(iniPar->GetValue()) );
			if(printTrue) truePar->SetValue( shiftAngle(truePar->GetValue()) );
		}

		// Is parameter a magnitude?
		bool isMag=0;
		if(outPar->GetName().find("mag")!=string::npos) isMag=1;
		// ... then make sure that it is positive
		if(isMag && !isFixed) {
			outPar->SetValue( std::fabs(outPar->GetValue()) );
			if(printInitial && iniPar)
				iniPar->SetValue( std::fabs(iniPar->GetValue()) );
			if(printTrue && truePar)
				truePar->SetValue( std::fabs(truePar->GetValue()) );
		}

		//Print parameter name
		*tableResult << o << outPar->GetName();

		//Print initial values
		if(printInitial)
			if( iniPar ) *tableResult << *iniPar;// |nr.| name| inital value|
			else *tableResult<<" ";

		//Print final value
		if( !isFixed ) *tableResult << *outPar;//final value
		else *tableResult<<" ";

		//Print true values
		if(printTrue)
			if( truePar ){
				*tableResult << *truePar;
				double pi = PhysConst::instance()->getConstValue("Pi");
				double pull = (truePar->GetValue()-outPar->GetValue() );
				// Shift pull by 2*pi if that reduces the deviation
				if(isAngle && !isFixed) {
					while( pull<0 && pull<-pi) pull+=2*pi;
					while( pull>0 && pull>pi) pull-=2*pi;
				}
				if(outPar->HasError()){
					if( errorType == ErrorType::ASYM && pull < 0)
						pull /= outPar->GetErrorLow();
					else if( errorType == ErrorType::ASYM && pull > 0)
						pull /= outPar->GetErrorHigh();
					else
						pull /= outPar->GetError();
				}
				if( !std::isnan(pull) )
					*tableResult << pull;
				else
					*tableResult << " ";
			} else *tableResult << " " << " ";
	}
	tableResult->footer();

	return;
}

void FitResult::printFitFractions(TableFormater* tab)
{
	BOOST_LOG_TRIVIAL(info) << " FitResult::printFitFractions() | "
			"Calculating fit fractions!";
	auto itrAmp = _ampVec.begin();
	for( ; itrAmp!=_ampVec.end(); ++itrAmp){
		printFitFractions(tab, (*itrAmp), nSetsFractionError);
	}
}

void FitResult::printFitFractions(TableFormater* fracTable,
		std::shared_ptr<Amplitude> amp, int nErrorSets)
{
	ParameterList ffList;
	amp->GetFitFractions(ffList);
	calcFractionError(ffList, amp, nErrorSets);
	double sum=0, sumErrorSq=0;

	fracTable->Reset();

	std::string ampName = amp->GetName();
	//print matrix
	fracTable->addColumn("Resonance: "+ampName,40);//add empty first column
	fracTable->addColumn("Fraction",15);//add empty first column
	fracTable->addColumn("Error",15);//add empty first column
	fracTable->header();
	for(unsigned int i=0;i<ffList.GetNDouble(); ++i){
		std::shared_ptr<DoubleParameter> tmpPar = ffList.GetDoubleParameter(i);
		std::string resName = tmpPar->GetName();

		//Remove amplitude name from string
		std::string::size_type strPos = resName.find(ampName);
		if (strPos != std::string::npos)
			resName.erase(strPos, ampName.length());

		*fracTable << resName
			<< tmpPar->GetValue()
			<< tmpPar->GetError(); //assume symmetric errors here
		sum += tmpPar->GetValue();
		sumErrorSq += tmpPar->GetError()*tmpPar->GetError();
	}
	fracTable->delim();
	*fracTable << "Total" << sum << sqrt(sumErrorSq) ;
	fracTable->footer();
	fractionList = ffList;
	sumFractions = sum;
	sumFractionsError = sqrt(sumErrorSq);

	return;
}

//void FitResult::calcFraction(ParameterList& parList, int nSets)
//{
//	if( !_ampVec.size() )
//		throw std::runtime_error("FitResult::calcFractions() | "
//				"No amplitude set, can't calculate fractions!");
//
//	if(parList.GetNDouble())
//		throw std::runtime_error("FitResult::calcFractions() | "
//				"ParameterList not empty!");
//
//	//	_amp->UpdateParameters(finalParameters); //update parameters in amplitude
//	double norm =-1;
//
//	//Start loop over amplitudes
//	auto ampItr = _ampVec.begin();
//	for( ; ampItr != _ampVec.end(); ++ampItr){
////		calcFraction(parList, (*ampItr));
//		(*ampItr)->GetFitFractions(parList);
//		if( parList.GetNDouble() )
//			calcFractionError(parList, (*ampItr), nSets);
//	}
//
//	return;
//}
} /* namespace ComPWA */


