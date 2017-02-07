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

#include <boost/timer.hpp>

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
#include "Core/FitResult.hpp"
//#include <boost/numeric/ublas/symmetric.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/log/core.hpp>
//#include <boost/log/trivial.hpp>
//#include <boost/timer.hpp>

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

using namespace boost::log;
using namespace ROOT::Minuit2;

double shiftAngle(double v)
{
	double originalVal = v;
	double val = originalVal;
	double pi = ComPWA::PhysConst::instance()->findConstant("Pi").value_;
	while(val> pi) val-=2*pi;
	while(val< -pi ) val+=2*pi;
	if(val!=originalVal)
		BOOST_LOG_TRIVIAL(info) << "shiftAngle() | Shifting parameter from "
		<<originalVal<< " to "<<val<<"!";
	return val;
}


MinuitIF::MinuitIF(std::shared_ptr<ControlParameter> esti, ParameterList& par) :
		_myFcn(esti, par), estimator(esti), enableHesse(true), enableMinos(true)
{

}

MinuitIF::~MinuitIF()
{

}

std::shared_ptr<FitResult> MinuitIF::exec(ParameterList& par)
{
	boost::timer time;
	par.RemoveDuplicates();

	ParameterList initialParList;
	initialParList.DeepCopy(par);

	MnUserParameters upar;
	int freePars = 0;
	for(unsigned int i=0; i<par.GetNDouble(); ++i){ //only doubles for minuit
		std::shared_ptr<DoubleParameter> actPat = par.GetDoubleParameter(i);
		//if no error is set or error set to 0 we use a default error,
		//otherwise minuit treads this parameter as fixed
		double error = actPat->GetError();
		if(error<=0) error = 0.01;
		if(!actPat->IsFixed() && actPat->GetName().find("phase") != actPat->GetName().npos )
			actPat->SetValue( shiftAngle(actPat->GetValue()) );

		if( actPat->UseBounds() ){
			upar.Add(actPat->GetName(), actPat->GetValue(), error,
					actPat->GetMinValue(), actPat->GetMaxValue());
		} else {
			upar.Add(actPat->GetName(), actPat->GetValue(),error);
		}

		if(!actPat->IsFixed())
			freePars++;
		if(actPat->IsFixed())
			upar.Fix(actPat->GetName());
	}

	BOOST_LOG_TRIVIAL(info) << "MinuitIF::exec() | Number of parameters (free): "
			<<par.GetNDouble()<<" ("<<freePars<<")";

	//use MnStrategy class to set all options for the fit
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
	BOOST_LOG_TRIVIAL(debug) << "Gradient number of steps: "
			<<strat.GradientNCycles();
	BOOST_LOG_TRIVIAL(debug) << "Gradient step tolerance: "
			<<strat.GradientStepTolerance();
	BOOST_LOG_TRIVIAL(debug) << "Gradient tolerance: "
			<<strat.GradientTolerance();
	BOOST_LOG_TRIVIAL(debug) << "Hesse number of steps: "
			<<strat.HessianNCycles();
	BOOST_LOG_TRIVIAL(debug) << "Hesse gradient number of steps: "
			<<strat.HessianGradientNCycles();
	BOOST_LOG_TRIVIAL(debug) << "Hesse step tolerance: "
			<<strat.HessianStepTolerance();
	BOOST_LOG_TRIVIAL(debug) << "Hesse G2 tolerance: "
			<<strat.HessianG2Tolerance();

	//MIGRAD
	MnMigrad migrad(_myFcn, upar, strat);
	double maxfcn = 0.0;
	double tolerance = 0.1;

	/* From the MINUIT2 Documentation:
	 * Minimize the function MnMigrad()(maxfcn, tolerance)
	 *    @param maxfcn : max number of function calls (if = 0) default is
	 *           used which is set to 200 + 100 * npar + 5 * npar**2
	 *    @param tolerance : value used for terminating iteration procedure.
	 *           For example, MIGRAD will stop iterating when edm (expected
	 *           distance from minimum) will be: edm < tolerance * 10**-3
	 *           Default value of tolerance used is 0.1
	 */
	BOOST_LOG_TRIVIAL(info) <<"MinuitIF::exec() | Starting migrad: "
			"maxCalls="<<maxfcn<<" tolerance="<<tolerance;

	FunctionMinimum minMin = migrad(maxfcn,tolerance);//(maxfcn,tolerance)

	BOOST_LOG_TRIVIAL(info) <<"MinuitIF::exec() | Migrad finished! "
			"Minimum is valid = "<<minMin.IsValid();

	//HESSE
	MnHesse hesse(strat);
	if( minMin.IsValid() && enableHesse ) {
		BOOST_LOG_TRIVIAL(info) <<"MinuitIF::exec() | Starting hesse";
		hesse(_myFcn,minMin);//function minimum minMin is updated by hesse
		BOOST_LOG_TRIVIAL(info) <<"MinuitIF::exec() | Hesse finished";
	} else
		BOOST_LOG_TRIVIAL(info) <<"MinuitIF::exec() | Migrad failed to "
				"find minimum! Skip hesse and minos!";

	BOOST_LOG_TRIVIAL(info) << "MinuitIF::exec() | Minimization finished! "
			"LH = "<<std::setprecision(10)<<minMin.Fval();

	//MINOS
	MnMinos minos(_myFcn,minMin,strat);

	//save minimzed values
	MnUserParameterState minState = minMin.UserState();

	//ParameterList can be changes by minos. We have to do a deep copy here.
	ParameterList finalParList;
	finalParList.DeepCopy(par);

	/* We directly write the central values of the fit result to check if the
	 * parameters change later on
	 */
	std::stringstream resultsOut;
	resultsOut<<"Central values of floating paramters:"<<std::endl;
	for(unsigned int i=0; i<finalParList.GetNDouble(); ++i){
		auto finalPar =	finalParList.GetDoubleParameter(i);
		if(finalPar->IsFixed()) continue;

		//central value
		double val = minState.Value(finalPar->GetName());

		//shift to [-pi;pi] if parameter is a phase
		if(finalPar->GetName().find("phase") != finalPar->GetName().npos)
			val =  shiftAngle(val);
		finalPar->SetValue(val);

		resultsOut<<finalPar->GetName()<<" "<<val<<std::endl;
		if( finalPar->GetErrorType()==ErrorType::ASYM  ){
			// Skip minos and fill symmetic errors
			if(!minMin.IsValid() || !enableMinos){
				BOOST_LOG_TRIVIAL(info) <<"MinuitIF::exec() | Skip Minos "
						"for parameter "<<i<< "...";
				finalPar->SetError(minState.Error(finalPar->GetName()));
				continue;
			}
			//asymmetric errors -> run minos
			BOOST_LOG_TRIVIAL(info) <<"MinuitIF::exec() | Run minos "
					"for parameter "<<i<< "...";
			MinosError err = minos.Minos(i);
			//lower = pair.first, upper= pair.second
			std::pair<double,double> assymErrors = err();
			finalPar->SetError( assymErrors.first, assymErrors.second );
		} else if(finalPar->GetErrorType()==ErrorType::SYM) {
			//symmetric errors -> migrad/hesse error
			finalPar->SetError(minState.Error(finalPar->GetName()));
		} else {
			throw std::runtime_error(
					"MinuitIF::exec() | Unknown error type of parameter: "
					+std::to_string(
							(long long int)finalPar->GetErrorType()
					)
			);
		}
	}
	BOOST_LOG_TRIVIAL(debug)<<"MinuitIF::exec() | "<<resultsOut.str();

	// Create fit result
	std::shared_ptr<FitResult> result(
			new MinuitResult(estimator, minMin)
	);
	result->setFinalParameters(finalParList);
	result->setInitialParameters(initialParList);
	result->setTime(time.elapsed());

	//update parameters in amplitude
	Amplitude::UpdateAmpParameterList(estimator->getAmplitudes(), finalParList);

	return result;
}

} /* namespace Minuit2 */
} /* namespace Optimizer */
} /* namespace ComPWA */
