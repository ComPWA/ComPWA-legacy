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
#include <vector>
#include <time.h>
#include <string>
#include <sstream>
#include <iostream>
#include <memory>
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MinosError.h"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/timer.hpp>
using namespace boost::log;
using namespace ROOT::Minuit2;

MinuitIF::MinuitIF(std::shared_ptr<ControlParameter> theData, ParameterList& par) : _myFcn(theData, par){
	//_myFcn = new MIMinuitFcn(theData);
}

MinuitIF::~MinuitIF(){
	//std::cout << "MinuitIF::~MinuitIF: I'll be back" << std::endl;
	//delete _myFcn;
}

//const double MinuitIF::exec(ParameterList& par){
std::shared_ptr<FitResult> MinuitIF::exec(ParameterList& par){
	boost::timer time;
	ParameterList initialParList(par);

	MnUserParameters upar;
	BOOST_LOG_TRIVIAL(debug) << "Parameters used: "<<par.GetNDouble();
	for(unsigned int i=0; i<par.GetNDouble(); ++i){ //only doubles for minuit

		//out.str("");
		//out << i;
		//s = out.str();

		//use as much information as possible: (just bounds but no error not supported by minuit)
		//try{
		std::shared_ptr<DoubleParameter> actPat = par.GetDoubleParameter(i);
		//}
		double error = 0.01;//if no error is set or error set to 0 we use a default error, otherwise minuit treads this parameter as fixed
		if( actPat->UseBounds() && actPat->HasError() ){
			if(actPat->GetError()!=0) error=*actPat->GetError();
			upar.Add(actPat->GetName(), actPat->GetValue(), error, actPat->GetMaxValue(), actPat->GetMinValue());
		}else if( actPat->HasError() ){
			if(actPat->GetError()!=0) error=*actPat->GetError();
			upar.Add(actPat->GetName(), actPat->GetValue(), error);
		}else
			upar.Add(actPat->GetName(), actPat->GetValue(),error);

		_myFcn.setNameID(i, actPat->GetName());

		if(actPat->IsFixed())
			upar.Fix(actPat->GetName());
	}

	//use MnStrategy class to set all options for the fit
	//	MnStrategy strat; //using default strategy = 1
	MinuitStrategy strat;//using default strategy = 1

	//read in xml configuration file for strategy settings
	const char* pPath = getenv("COMPWA_DIR");
	std::string path = std::string(pPath);
	std::ifstream ifs(path+"/Optimizer/Minuit2/MinuitStrategy.xml");
	boost::archive::xml_iarchive ia(ifs,boost::archive::no_header);
	ia >> BOOST_SERIALIZATION_NVP(strat);
	strat.init();//update parameters of MnStrategy mother class (IMPORTANT!)
	ifs.close();
	//write strategy settings
	//	std::ofstream ofs(path+"Optimizer/Minuit2/test.xml");
	//	boost::archive::xml_oarchive oa(ofs,boost::archive::no_header);
	//	oa << BOOST_SERIALIZATION_NVP(strat);
	//	ofs.close();
	BOOST_LOG_TRIVIAL(debug) << "Minuit strategy parameters: (from "<<
			path+"Optimizer/Minuit2/MinuitStrategy.xml"<<")";
	BOOST_LOG_TRIVIAL(debug) << "Gradient number of steps: "<<strat.GradientNCycles();
	BOOST_LOG_TRIVIAL(debug) << "Gradient step tolerance: "<<strat.GradientStepTolerance();
	BOOST_LOG_TRIVIAL(debug) << "Gradient tolerance: "<<strat.GradientTolerance();
	BOOST_LOG_TRIVIAL(debug) << "Hesse number of steps: "<<strat.HessianNCycles();
	BOOST_LOG_TRIVIAL(debug) << "Hesse gradient number of steps: "<<strat.HessianGradientNCycles();
	BOOST_LOG_TRIVIAL(debug) << "Hesse step tolerance: "<<strat.HessianStepTolerance();
	BOOST_LOG_TRIVIAL(debug) << "Hesse G2 tolerance: "<<strat.HessianG2Tolerance();

	//MIGRAD
	BOOST_LOG_TRIVIAL(info) <<"MinuitIF: starting migrad ";
	MnMigrad migrad(_myFcn, upar, strat);
	//	FunctionMinimum minMin = migrad(100,0.001);//(maxfcn,tolerance)
	FunctionMinimum minMin = migrad();
	BOOST_LOG_TRIVIAL(info) <<"MinuitIF: migrad finished";

	//HESSE
	BOOST_LOG_TRIVIAL(info) <<"MinuitIF: starting hesse";
	MnHesse hesse(strat);
	hesse(_myFcn,minMin);//function minimum minMin is updated by hesse
	BOOST_LOG_TRIVIAL(info) <<"MinuitIF: hesse finished";

	//MINOS
	MnMinos minos(_myFcn,minMin,strat);

	//we copy parameters here because minos can still change the parameterList par
	ParameterList finalParList(par);
	//save minimzed values
	MnUserParameterState minState = minMin.UserState();
	for(unsigned int i=0; i<finalParList.GetNDouble(); ++i){
		std::shared_ptr<DoubleParameter> finalPar = finalParList.GetDoubleParameter(i);
		if(!finalPar->IsFixed()){
			finalPar->SetValue(minState.Value(finalPar->GetName()));
			if(finalPar->GetErrorType()==ErrorType::ASYM){ //asymmetric errors -> run minos
				BOOST_LOG_TRIVIAL(info) <<"MinuitIF: minos for parameter "<<i<< "...";
				MinosError err = minos.Minos(i);
				std::pair<double,double> assymErrors = err();//lower = pair.first, upper= pair.second
				finalPar->SetError( std::shared_ptr<ParError<double>>(new AsymError<double>(assymErrors)) );
			} else if(finalPar->GetErrorType()==ErrorType::SYM) {//symmetric errors -> migrad error
				finalPar->SetError(std::shared_ptr<ParError<double>>(
						new SymError<double>(minState.Error(finalPar->GetName()))));
			} else {
				BOOST_LOG_TRIVIAL(error)<< "MinuitIF: requesting error type "<<finalPar->GetErrorType()<<". No idea what to do here!";
				exit(1);
			}
		}
	}
//	std::cout<<"34234748 ";
//	for(unsigned int i=0; i<finalParList.GetNDouble(); ++i){
//		std::shared_ptr<DoubleParameter> finalPar = finalParList.GetDoubleParameter(i);
//		if(!finalPar->IsFixed()){
//			std::cout<<finalPar->GetName()
//							<<" "<<minState.Value(finalPar->GetName())
//							<<" "<<finalPar->GetValue()
//							<<" "<<minState.Error(finalPar->GetName())
//							<<" "<<finalPar->GetError()->GetErrorLow()
//							<<" "<<finalPar->GetError()->GetErrorHigh()<<" ";
//		}
//	}
//	std::cout<<std::endl;

	std::shared_ptr<FitResult> result(new MinuitResult(minMin));
	result->setInitialParameters(initialParList);
	result->setFinalParameters(finalParList);
	result->setTime(time.elapsed());

	return result;
}

