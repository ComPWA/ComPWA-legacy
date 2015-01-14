/*
 * NonResonant.hpp
 *
 *  Created on: Jan 13, 2015
 *      Author: weidenka
 */

#ifndef NONRESONANT_HPP_
#define NONRESONANT_HPP_

#include "Core/DataPoint.hpp"
#include "Core/Parameter.hpp"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpKinematics.hpp"

class NonResonant : public AmpAbsDynamicalFunction {
public:
	NonResonant(std::string name);
	virtual void initialise() { };
	//! value at \param point
	virtual std::complex<double> evaluate(dataPoint& point) { return _norm*evaluateAmp(point); };
	//! value of dynamical amplitude at \param point
	virtual std::complex<double> evaluateAmp(dataPoint& point) { return 1;} ;
	//! value of WignerD amplitude at \param point
	virtual double evaluateWignerD(dataPoint& point) { return 1;} ;

	//! Get resonance spin
	virtual double getSpin() { return 0;};

	static std::complex<double> dynamicalFunction();
protected:

};

class NonResonantStrategy : public Strategy {
public:
	NonResonantStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}

	virtual const std::string to_str() const {
		return ("NonResonant "+name);
	}

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
		if( checkType != out->type() ) {
			throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
			return false;
		}
		//MultiDim output, must have multidim Paras in input
		if(checkType == ParType::MCOMPLEX){
			if(paras.GetNMultiDouble()){
				unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();
				std::vector<std::complex<double> > results(nElements, std::complex<double>(0.));
				for(unsigned int ele=0; ele<nElements; ele++){
					results[ele] = NonResonant::dynamicalFunction();
				}
				out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));
				return true;
			}else{ //end multidim para treatment
				throw(WrongParType("Input MultiDoubles missing in BW strat of "+name));
				return false;
			}
		}

		std::complex<double> result = NonResonant::dynamicalFunction();
		out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
		return true;
	}

protected:
	std::string name;
};




#endif /* NONRESONANT_HPP_ */
