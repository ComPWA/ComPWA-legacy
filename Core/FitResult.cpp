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
	double pi = PhysConst::Instance().findConstant("Pi").value_;
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
		if(printInitial){
			try{
				iniPar = initialParameters.GetDoubleParameter(o);
			} catch (BadParameter& bad){
				BOOST_LOG_TRIVIAL(error) << "FitResult::printFitParameters() | "
						"can't access parameter of initial parameter list!";
				throw;
			}
		}
		try{
			outPar = finalParameters.GetDoubleParameter(o);
		} catch (BadParameter& bad){
			BOOST_LOG_TRIVIAL(error) << "FitResult::printFitParameters() | "
					"can't access parameter of final parameter list!";
			throw;
		}
		ErrorType errorType = outPar->GetErrorType();
		bool isFixed = outPar->IsFixed();
		bool isAngle=0;
		if(outPar->GetName().find("phase")!=string::npos) isAngle=1;//is our Parameter an angle?
		if(isAngle && !isFixed) {
			outPar->SetValue( shiftAngle(outPar->GetValue()) ); //shift angle to the interval [-pi;pi]
		}

		*tableResult << o << outPar->GetName();
		if(printInitial) *tableResult << *iniPar;// |nr.| name| inital value|
		if(isFixed) *tableResult<<"FIXED";
		else
			*tableResult << *outPar;//final value
		if(printTrue){
			try{
				truePar = trueParameters.GetDoubleParameter(outPar->GetName());
			} catch (BadParameter& bad){
				BOOST_LOG_TRIVIAL(error) << "FitResult::printFitParameters() | "
						"can't access parameter of true parameter list!";
				*tableResult << "not found"<< " - ";
				continue;
			}
			*tableResult << *truePar;
			double pi = PhysConst::Instance().findConstant("Pi").value_;
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
