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
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
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
	//std::string s;
	//std::stringstream out;

	unsigned int startTime = clock();
	ParameterList initialPar(par);

	MnUserParameters upar;
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
			if(actPat->GetError()!=0) error=actPat->GetError();
			upar.Add(actPat->GetName(), actPat->GetValue(), error, actPat->GetMaxValue(), actPat->GetMinValue());
		}else if( actPat->HasError() ){
			if(actPat->GetError()!=0) error=actPat->GetError();
			upar.Add(actPat->GetName(), actPat->GetValue(), error);
		}else
			upar.Add(actPat->GetName(), actPat->GetValue(),error);

		_myFcn.setNameID(i, actPat->GetName());

		if(actPat->IsFixed())
			upar.Fix(actPat->GetName());
	}
	//use MnStrategy class to set all options for the fit
	MnStrategy strat; //using default strategy

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
	BOOST_LOG_TRIVIAL(info) <<"MinuitIF: starting minos";
	MnMinos minos(_myFcn,minMin,strat);
	BOOST_LOG_TRIVIAL(info) <<"MinuitIF: minos for parameter 1 ...";
//	minos(1); //first parameter
	BOOST_LOG_TRIVIAL(info) <<"MinuitIF: minos finished";

	//save minimzed values
	MnUserParameterState minState = minMin.UserState();
	for(unsigned int i=0; i<par.GetNDouble(); ++i){
		std::shared_ptr<DoubleParameter> actPat = par.GetDoubleParameter(i);
		if(!actPat->IsFixed()){
			actPat->SetValue(minState.Value(actPat->GetName()));
			actPat->SetError(minState.Error(actPat->GetName()));
		}
	}
	std::shared_ptr<FitResult> result(new MinuitResult(minMin));
	result->setInitialParameters(initialPar);
	result->setFinalParameters(par);
	result->setTime((clock()-startTime)/CLOCKS_PER_SEC);

	return result;
}

