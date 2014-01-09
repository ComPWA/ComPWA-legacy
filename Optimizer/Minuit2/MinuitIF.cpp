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
	std::shared_ptr<FitResult> result(new FitResult());
	result->initialParameters=ParameterList(par);
	result->initialLH=-999; //how can i get to initial LH?
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
		}else{
			//double tmpVal = actPat->GetValue();
			//if(tmpVal==0 && !actPat->IsFixed()) tmpVal=0.001;
			upar.Add(actPat->GetName(), actPat->GetValue(),error);
		}

		_myFcn.setNameID(i, actPat->GetName());

		if(actPat->IsFixed())
			upar.Fix(actPat->GetName());
	}

	MnMigrad migrad(_myFcn, upar);
	BOOST_LOG_TRIVIAL(info) <<"MinuitIF: starting migrad ";
	//for(unsigned int i=0; i<par.GetNDouble(); i++)
	//   std::cout << upar.Parameter(i).Value() << " " << upar.Parameter(i).IsFixed() << std::endl;
	FunctionMinimum minMin = migrad(100,0.001);//TODO

	//if(!minMin.IsValid()) {
	//try with higher strategy
	//    std::cout <<"FM is invalid, try with strategy = 2."<< std::endl;
	//   MnMigrad migrad2(_myFcn, minMin.UserState(), MnStrategy(2));
	//   minMin = migrad2(10,0.1);//TODO
	// }

	//save minimzed values
	for(unsigned int i=0; i<par.GetNDouble(); ++i){
		//out.str("");
		//out << i;
		// s = out.str();
		std::shared_ptr<DoubleParameter> actPat = par.GetDoubleParameter(i);
		if(!actPat->IsFixed()){
			actPat->SetValue(minMin.UserState().Value(actPat->GetName()));
			actPat->SetError(minMin.UserState().Error(actPat->GetName()));
		}
	}
	std::vector<double> minuitCovM = minMin.UserState().Covariance().Data();//Covariance matrix is empty !?
//	std::cout<<minuitCovM.size()<<" "<<minMin.UserState().Covariance().size()<<std::endl;
	unsigned int matrixSize = par.GetNParameter();
	using namespace boost::numeric::ublas;
	boost::numeric::ublas::symmetric_matrix<double,boost::numeric::ublas::upper> covMatrix(matrixSize,matrixSize);
	if(minuitCovM.size()==matrixSize*(matrixSize+1)/2){
		for (unsigned i = 0; i < covMatrix.size1 (); ++ i)
			for (unsigned j = i; j < covMatrix.size2 (); ++ j)
				covMatrix (i, j) = minuitCovM[3*i+j];
		std::cout << covMatrix << std::endl;
	} else BOOST_LOG_TRIVIAL(info)<<"MinuitIF: no valid correlation matrix available!";
	result->cov=covMatrix;
	result->finalParameters=ParameterList(par);
	result->finalLH = minMin.Fval();
	result->edm= minMin.Edm();
	result->time = (clock()-startTime)/CLOCKS_PER_SEC;

	//	return minMin.Fval();
	return result;
}

