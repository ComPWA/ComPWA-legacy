/*
 * NonResonant.hpp
 *
 *  Created on: Jan 13, 2015
 *      Author: weidenka
 */

#ifndef NONRESONANT_HPP_
#define NONRESONANT_HPP_

#include <boost/property_tree/ptree.hpp>

#include "Core/DataPoint.hpp"
#include "Core/Parameter.hpp"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"

using boost::property_tree::ptree;

class NonResonant : public AmpAbsDynamicalFunction {
public:
	NonResonant(const char *name,
			std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
			int nCalls=30000, normStyle nS=normStyle::one) ;

	virtual void initialise() { };
	//! value of dynamical amplitude at \param point
	virtual std::complex<double> evaluateAmp(dataPoint& point) { return dynamicalFunction();}
	//! value of WignerD amplitude at \param point
	virtual double evaluateWignerD(dataPoint& point) { return 1;} ;
	//!Integral
	virtual double integral() { return 1/sqrt(Kinematics::instance()->getPhspVolume()); }

	static std::complex<double> dynamicalFunction();

	virtual std::shared_ptr<FunctionTree> setupTree(
			allMasses& theMasses,allMasses& toyPhspSample,std::string suffix, ParameterList& params);
protected:

};

class basicConf
{
public:
	virtual ~basicConf(){};
	basicConf(){}
	basicConf(const boost::property_tree::ptree &pt_){
		m_enable = pt_.get<bool>("enable");
		m_name= pt_.get<std::string>("name");
		m_strength= pt_.get<double>("strength");
		m_strength_fix = pt_.get<bool>("strength_fix");
		m_strength_min= pt_.get<double>("strength_min");
		m_strength_max= pt_.get<double>("strength_max");
		m_phase= pt_.get<double>("phase");
		m_phase_fix= pt_.get<bool>("phase_fix");
		m_phase_min= pt_.get<double>("phase_min");
		m_phase_max= pt_.get<double>("phase_max");
		try{
			m_ffType= pt_.get<int>("FormFactorType");
		} catch (...) {
			m_ffType = 1;
		}
	}
	virtual void put(boost::property_tree::ptree &pt_){
		pt_.put("enable", m_enable);
		pt_.put("name", m_name);
		pt_.put("strength", m_strength);
		pt_.put("strength_fix", m_strength_fix);
		pt_.put("strength_min", m_strength_min);
		pt_.put("strength_max", m_strength_max);
		pt_.put("phase", m_phase);
		pt_.put("phase_fix", m_phase_fix);
		pt_.put("phase_min", m_phase_min);
		pt_.put("phase_max", m_phase_max);
		pt_.put("FormFactorType", m_ffType);
	}
	virtual void update(ParameterList par){
		try{// only update parameters if they are found in list
			m_strength= par.GetDoubleParameter("mag_"+m_name)->GetValue();
			m_strength_fix= par.GetDoubleParameter("mag_"+m_name)->IsFixed();
			m_strength_min= par.GetDoubleParameter("mag_"+m_name)->GetMinValue();
			m_strength_max= par.GetDoubleParameter("mag_"+m_name)->GetMaxValue();
			m_phase = par.GetDoubleParameter("phase_"+m_name)->GetValue();
			m_phase_fix= par.GetDoubleParameter("phase_"+m_name)->IsFixed();
			m_phase_min= par.GetDoubleParameter("phase_"+m_name)->GetMinValue();
			m_phase_max= par.GetDoubleParameter("phase_"+m_name)->GetMaxValue();
			m_ffType= par.GetDoubleParameter("ffType_"+m_name)->GetValue();
		} catch (BadParameter b) {
			//BOOST_LOG_TRIVIAL(error)<<"basicConf::update() | one or more parameters haven't been found!";
		}
	}
	//protected:
	bool m_enable;
	std::string m_name;
	double m_strength;
	bool m_strength_fix;
	double m_strength_min;
	double m_strength_max;
	double m_phase;
	bool m_phase_fix;
	double m_phase_min;
	double m_phase_max;
	int m_ffType;
};




#endif /* NONRESONANT_HPP_ */
