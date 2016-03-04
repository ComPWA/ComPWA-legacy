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

using boost::property_tree::ptree;

class NonResonant : public AmpAbsDynamicalFunction
{
public:

	NonResonant( normStyle nS=normStyle::one, int calls=30000 ) :
		AmpAbsDynamicalFunction( nS, calls ) { };

	NonResonant(const char *name,
			std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
			int nCalls=30000, normStyle nS=normStyle::one) ;

	virtual void initialise() { };
	//! Configure resonance from ptree
	virtual void Configure(boost::property_tree::ptree::value_type const& v,
			ParameterList& list);

	virtual void Save(boost::property_tree::ptree &pt);
	//! Get resonance width
	virtual double GetWidth() const { return 0; }
	//! value of dynamical amplitude at \param point
	virtual std::complex<double> EvaluateAmp(dataPoint& point) {
		return dynamicalFunction();
	}
	//! value of WignerD amplitude at \param point
	virtual double EvaluateWignerD(dataPoint& point) { return 1;} ;
	//!Integral
	virtual double GetIntegral() {
		return 1/sqrt(Kinematics::instance()->getPhspVolume());
	}

	static std::complex<double> dynamicalFunction();

	virtual std::shared_ptr<FunctionTree> SetupTree(
			allMasses& theMasses,allMasses& toyPhspSample,std::string suffix);
};
#endif /* NONRESONANT_HPP_ */
