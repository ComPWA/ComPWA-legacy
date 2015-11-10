/*
 * FitResult.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: weidenka
 */
#include "Core/FitResult.hpp"

void FitResult::writeText(std::string filename){
	std::ofstream myfile;
	myfile.open(filename);
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
	double pi = PhysConst::instance()->getConstValue("Pi");
	while(val> pi) val-=2*pi;
	while(val< -pi ) val+=2*pi;
	if(val!=originalVal)
		BOOST_LOG_TRIVIAL(info) << "shiftAngle(): shifting parameter from "<<originalVal<< " to "<<val<<"!";
	return val;
}

void FitResult::genSimpleOutput(std::ostream& out){
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> outPar = finalParameters.GetDoubleParameter(o);
		out<<outPar->GetValue()<<" "<<outPar->GetError()<<" ";
	}
	out<<"\n";

	return;
}
void FitResult::setFinalParameters(ParameterList finPars){
	finalParameters=finPars;
	_amp->setParameterList(finalParameters); //update parameters in amplitude
}

void FitResult::printFitParameters(TableFormater* tableResult){
	bool printTrue=0, printInitial=0;
	if(trueParameters.GetNParameter()) printTrue=1;
	if(initialParameters.GetNParameter()) printInitial=1;
	unsigned int parErrorWidth = 22;
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++)
		if(finalParameters.GetDoubleParameter(o)->GetErrorType()==ErrorType::ASYM) parErrorWidth=33;

	tableResult->addColumn("Nr");
	tableResult->addColumn("Name",15);
	if(printInitial) tableResult->addColumn("Initial Value",parErrorWidth);
	tableResult->addColumn("Final Value",parErrorWidth);
	if(printTrue) tableResult->addColumn("True Value",13);
	if(printTrue) tableResult->addColumn("Deviation",13);
	tableResult->header();

	std::shared_ptr<DoubleParameter> iniPar, outPar, truePar;
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
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
				BOOST_LOG_TRIVIAL(error) << "FitResult::printFitParameters() | "
						"can't access parameter '"<<outPar->GetName()<<"' of initial parameter list!";
				throw;
			}
		}
		if(printTrue){
			try{
				truePar = trueParameters.GetDoubleParameter(outPar->GetName());
			} catch (BadParameter& bad){
				BOOST_LOG_TRIVIAL(error) << "FitResult::printFitParameters() | "
						"can't access parameter '"<<outPar->GetName()<<"' of true parameter list!";
				*tableResult << "not found"<< " - ";
				throw;
			}
		}

		ErrorType errorType = outPar->GetErrorType();
		bool isFixed = outPar->IsFixed();
		bool isAngle=0, isMag=0;
		if(outPar->GetName().find("phase")!=string::npos) isAngle=1;//is our Parameter an angle?
		if(outPar->GetName().find("mag")!=string::npos) isMag=1;//is our Parameter an angle?
		if(isAngle && !isFixed) {
			outPar->SetValue( shiftAngle(outPar->GetValue()) ); //shift angle to the interval [-pi;pi]
			if(printInitial) iniPar->SetValue( shiftAngle(iniPar->GetValue()) );
			if(printTrue) truePar->SetValue( shiftAngle(truePar->GetValue()) );
		}
		if(isMag && !isFixed) {
			outPar->SetValue( std::abs(outPar->GetValue()) ); //abs value of parameter is magnitude
			if(printInitial) iniPar->SetValue( std::abs(iniPar->GetValue()) );
			if(printTrue) truePar->SetValue( std::abs(truePar->GetValue()) );
		}

		*tableResult << o << outPar->GetName();
		if(printInitial) *tableResult << *iniPar;// |nr.| name| inital value|
		if(isFixed) *tableResult<<"FIXED";
		else
			*tableResult << *outPar;//final value
		if(printTrue){
			*tableResult << *truePar;
			double pi = PhysConst::instance()->getConstValue("Pi");
			double pull = (truePar->GetValue()-outPar->GetValue() );
			if(isAngle && !isFixed) { //shift pull by 2*pi if that reduces the deviation
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
			*tableResult << pull;
		}
	}
	tableResult->footer();

	return;
}

void FitResult::printFitFractions(TableFormater* fracTable){
	BOOST_LOG_TRIVIAL(info) << "Calculating fit fractions...";
	calcFraction();
	double sum, sumErrorSq;

	//print matrix
	fracTable->addColumn("Resonance",15);//add empty first column
	fracTable->addColumn("Fraction",15);//add empty first column
	fracTable->addColumn("Error",15);//add empty first column
	fracTable->addColumn("Significance",15);//add empty first column
	fracTable->header();
	for(unsigned int i=0;i<fractionList.GetNDouble(); ++i){
		std::shared_ptr<DoubleParameter> tmpPar = fractionList.GetDoubleParameter(i);
		*fracTable << tmpPar->GetName()
				<< tmpPar->GetValue()
				<< tmpPar->GetError() //assume symmetric errors here
				<< std::abs(tmpPar->GetValue()/tmpPar->GetError());
		sum += tmpPar->GetValue();
		sumErrorSq += tmpPar->GetError()*tmpPar->GetError();
	}
	fracTable->delim();
	*fracTable << "Total" << sum << sqrt(sumErrorSq) << " ";
	fracTable->footer();

	return;
}

void FitResult::calcFraction() {
	if(!fractionList.GetNDouble()) {
		calcFraction(fractionList);
		calcFractionError();
	} else
		BOOST_LOG_TRIVIAL(warning) << "FitResult::calcFractions() fractions already calculated. Skip!";
}

void FitResult::calcFraction(ParameterList& parList){
	if(!_amp)
		throw std::runtime_error("FitResult::calcFractions() | no amplitude set, can't calculate fractions!");
	if(parList.GetNDouble())
		throw std::runtime_error("FitResult::calcFractions() | ParameterList not empty!");

	double norm =-1;
	ParameterList cList;
	_amp->copyParameterList(cList);

	//in case of unbinned efficiency correction to tree does not provide an integral w/o efficiency correction
	norm = _amp->integral();
	if(norm<0)
		throw std::runtime_error("FitResult::calcFraction() normalization can't be calculated");
	BOOST_LOG_TRIVIAL(debug)<<"FitResult::calcFraction() norm="<<norm;
	int nRes=_amp->GetNumberOfResonances();
	for(unsigned int i=0;i<nRes; i++){ //fill matrix
		double resInt= _amp->GetIntegral(i);
		std::string resName = _amp->GetNameOfResonance(i);
		std::shared_ptr<DoubleParameter> magPar = cList.GetDoubleParameter("mag_"+resName);
		double mag = std::abs(magPar->GetValue()); //value of magnitude
		double magError = 0;
		if(magPar->HasError()) magError = magPar->GetError(); //error of magnitude
		parList.AddParameter(std::shared_ptr<DoubleParameter>(
				new DoubleParameter(resName+"_FF", mag*mag*resInt/norm, 2*mag*resInt/norm * magError)) );
	}
	return;
}
