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
#include <string>
#include <sstream>
#include <iostream>

#include "Physics/BreitWigner/BreitWigner.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"

using namespace std;

BreitWigner::BreitWigner(const double min, const double max):min_(min),max_(max) {
	result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BreitWignerResult")));
	//_myFcn = new MIMinuitFcn(theData);
}

BreitWigner::~BreitWigner()
{
	//delete _myFcn;
}

const double BreitWigner::integral(ParameterList& par){
	double integral = 0;
	unsigned int nSteps = 1000000;
	double step = (max_-min_)/(double)nSteps;

	//TODO: try exceptions, parameter find via name?
	integral += step*BreitWignerValue(min_, par.GetDoubleParameter(0)->GetValue(), par.GetDoubleParameter(1)->GetValue())/2.;
	for(unsigned int k=1; k<nSteps; k++){
		integral += step*BreitWignerValue((min_+k*step), par.GetDoubleParameter(0)->GetValue(), par.GetDoubleParameter(1)->GetValue());
	}
	integral += step*BreitWignerValue(max_, par.GetDoubleParameter(0)->GetValue(), par.GetDoubleParameter(1)->GetValue())/2.;

	return integral;
}

const double BreitWigner::volume(){
	return 0;
}

const double BreitWigner::drawInt(double* x, double *p){
	return p[2]*BreitWignerValue(x[0], p[0], p[1]);
}

const ParameterList& BreitWigner::intensity(double x, double M, double T){
	//ParameterList result;
	//result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BreitWignerResult",BreitWignerValue(x, M, T))));
	result.SetParameterValue(0,BreitWignerValue(x, M, T));
	return result;
}

const ParameterList& BreitWigner::intensity(const dataPoint& point){
	//ParameterList result;
	double val = BreitWignerValue(point.getVal(0),params.GetDoubleParameter(0)->GetValue(), params.GetDoubleParameter(1)->GetValue());
	//result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BreitWignerResult",val)));
	result.SetParameterValue(0,val);
	return result;
}
const ParameterList& BreitWigner::intensity(const dataPoint& point, ParameterList& par){
	setParameterList(par);
	return intensity(point);
}
const ParameterList& BreitWigner::intensity(std::vector<double> x, ParameterList& par){
	//assmue that x[0]=mass23sq;
	dataPoint dataP; dataP.setVal(0,x[0]); dataP.setVal(1,x[1]);
	return intensity(dataP,par);
}

bool BreitWigner::copyParameterList(ParameterList& outPar){
	if(outPar.GetNParameter())
		return false; //already filled ,TODO: exception?

	outPar.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BreitWignerPosition",1.5, 0.5, 2.5, 0.1)));
	outPar.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BreitWignerWidth",0.3, 0.1, 0.5, 0.05)));

	return true;
}

const double BreitWigner::BreitWignerValue(double x, double M, double T){
	double denom=(x*x-M*M)*(x*x-M*M)+M*M*T*T;

	return 1./denom;
}
void BreitWigner::setParameterList(ParameterList& par){
	//parameters varied by Minimization algorithm
	params = par;
	return;
}
