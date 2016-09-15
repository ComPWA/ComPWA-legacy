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

namespace ComPWA {
namespace Physics {
namespace BreitWigner {

BreitWigner::BreitWigner(const double min, const double max):min_(min),max_(max)
{
	result.AddParameter(
			std::shared_ptr<DoubleParameter>(
					new DoubleParameter("BreitWignerResult")
			)
	);
}

BreitWigner::~BreitWigner()
{
	//delete _myFcn;
}

const double BreitWigner::GetIntegral()
{
	double integral = 0;
	unsigned int nSteps = 1000000;
	double step = (max_-min_)/(double)nSteps;

	//TODO: try exceptions, parameter find via name?
	integral += 0.5*step*BreitWignerValue(
			min_,
			params.GetDoubleParameter(0)->GetValue(),
			params.GetDoubleParameter(1)->GetValue()
	);

	for(unsigned int k=1; k<nSteps; k++){
		integral += step*BreitWignerValue(
				(min_+k*step),
				params.GetDoubleParameter(0)->GetValue(),
				params.GetDoubleParameter(1)->GetValue()
		);
	}

	integral += 0.5*step*BreitWignerValue(
			max_,
			params.GetDoubleParameter(0)->GetValue(),
			params.GetDoubleParameter(1)->GetValue()
	);

	return integral;
}

const double BreitWigner::volume()
{
	return 0;
}

const double BreitWigner::drawInt(double* x, double *p)
{
	return p[2]*BreitWignerValue(x[0], p[0], p[1]);
}

const ParameterList& BreitWigner::intensity(double x, double M, double T)
{
	//ParameterList result;
	//result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BreitWignerResult",BreitWignerValue(x, M, T))));
	result.SetParameterValue(0,BreitWignerValue(x, M, T));
	return result;
}

const ParameterList& BreitWigner::intensity(const dataPoint& point)
{
	//ParameterList result;
	double val = BreitWignerValue(
			point.getVal(0),
			params.GetDoubleParameter(0)->GetValue(),
			params.GetDoubleParameter(1)->GetValue()
			);
	result.SetParameterValue(0,val);
	return result;
}

const ParameterList& BreitWigner::intensity(std::vector<double> x)
{
	dataPoint point;
	try{
		Kinematics::instance()->FillDataPoint(1,0,x[1],x[0],point);
	} catch (BeyondPhsp& ex){
		result.SetParameterValue(0,0);
		return result;
	}
	return intensity(point);
}

const double BreitWigner::BreitWignerValue(double x, double M, double T){
	double denom=(x*x-M*M)*(x*x-M*M)+M*M*T*T;

	return 1./denom;
}

//void BreitWigner::setParameterList(const ParameterList& par){
//	//parameters varied by Minimization algorithm
//	params = par;
//	return;
//}

} /* namespace BreitWigner */
} /* namespace Physics */
} /* namespace ComPWA */
