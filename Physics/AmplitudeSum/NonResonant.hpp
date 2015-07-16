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
#include "Physics/AmplitudeSum/AmpKinematics.hpp"

using boost::property_tree::ptree;

class NonResonant : public AmpAbsDynamicalFunction {
public:
	NonResonant(std::string name);
	virtual void initialise() { };
	//! value at \param point
	virtual std::complex<double> evaluate(dataPoint& point) {
		if(GetNormalization()<0) return evaluateAmp(point); //normalization is disabled
		return (1/Kinematics::instance()->getPhspVolume())*evaluateAmp(point);
	}
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
	NonResonantStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){
	}

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
	}
	virtual void update(ParameterList par){
		try{// only update parameters if they are found in list
			m_strength= par.GetDoubleParameter("mag_"+m_name)->GetValue();
			m_phase = par.GetDoubleParameter("phase_"+m_name)->GetValue();
		} catch (BadParameter b) { }
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
};




#endif /* NONRESONANT_HPP_ */
