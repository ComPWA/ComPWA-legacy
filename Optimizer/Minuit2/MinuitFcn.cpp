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
#include <cassert>
#include <memory>
#include <iostream>
#include <cmath>

#include <boost/chrono.hpp>

#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"
#include "Optimizer/Minuit2/MinuitFcn.hpp"
#include "Optimizer/ControlParameter.hpp"
using namespace boost::log;
using namespace ROOT::Minuit2;

MinuitFcn::MinuitFcn(std::shared_ptr<ControlParameter> myData, ParameterList& parList) :
		  _myDataPtr(myData), _parList(parList){
	if (0==_myDataPtr) {
		//Alert << "Data pointer is 0 !!!!" << endmsg;
		std::cout << "Data pointer is 0 !!!!" << std::endl; //TODO exception
		exit(1);
	}
}

MinuitFcn::~MinuitFcn(){
	//std::cout << "~MinuitFcn: I'll be back" << std::endl;
}

double MinuitFcn::operator()(const std::vector<double>& x) const{
	//ParameterList par;
	std::ostringstream paramOut;
	for(unsigned int i=0; i<x.size(); i++){
		//par.AddParameter(DoubleParameter(std::string("tmpPar"+i),x[i]));
		//_parList.SetParameterValue(i,x[i]);
		std::shared_ptr<DoubleParameter> actPat = _parList.GetDoubleParameter(_parNames.at(i));
		if(!actPat->IsFixed())
			if(x[i]==x[i]){
				actPat->SetValue(x[i]);
				paramOut << x[i] << " ";//print only free parameters
//				std::cout<<_parNames.at(i)<<" "<<x[i]<<std::endl;
			}
	}
	boost::chrono::system_clock::time_point start = boost::chrono::system_clock::now();
	double result=_myDataPtr->controlParameter(_parList);
	boost::chrono::duration<double> sec = boost::chrono::system_clock::now() - start;
	BOOST_LOG_TRIVIAL(debug) << std::setprecision(10)
	<< "MinuitFcn: -log(L) = "<< result
	<< std::setprecision(4)
	<<" Time: "<<sec.count()<<"s"
	<<" nCalls: "<<_myDataPtr->nCalls()<<std::endl
	<<" Parameters: "<<paramOut.str();

	return result;
}

double MinuitFcn::Up() const{
	return 0.5; //TODO: Setter, LH 0.5, Chi2 1.
}



