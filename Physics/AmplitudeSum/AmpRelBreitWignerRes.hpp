//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef AMP_REL_BREIT_WIGNER_RES
#define AMP_REL_BREIT_WIGNER_RES

#include <vector>

#include "Core/Functions.hpp"
#include "Core/Exceptions.hpp"
#include "Core/FunctionTree.hpp"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Physics/AmplitudeSum/NonResonant.hpp"

class AmpRelBreitWignerRes : public AmpAbsDynamicalFunction {
public:

	AmpRelBreitWignerRes(const char *name,
			std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
			std::shared_ptr<DoubleParameter> mass, int subSys, Spin spin, Spin m, Spin n,
			std::shared_ptr<DoubleParameter> width,
			std::shared_ptr<DoubleParameter> mesonRadius,
			std::shared_ptr<DoubleParameter> motherRadius,
			formFactorType type = formFactorType::BlattWeisskopf,
			int nCalls=30000, normStyle nS=normStyle::one) ;

	virtual ~AmpRelBreitWignerRes();
	/** Breit-Wigner function
	 *
	 * The dynamical function implemented here is taken from PDG2014 (Eq.47-22) for
	 * the one channel case.
	 * @J is angular momentum between A&B. In case A&B have spin 0 this is the spin for the resonace.
	 *
	 * @param mSq Invariant mass
	 * @param mR Resonance mass
	 * @param ma Mass particle A
	 * @param mb Mass particle B
	 * @param width  Width of resonance
	 * @param J Angular momentum between A&B
	 * @param mesonRadius Scale of interaction range
	 * @return
	 */
	static std::complex<double> dynamicalFunction(double mSq, double mR, double ma, double mb,
			double width, unsigned int J, double mesonRadius,
			formFactorType ffType=formFactorType::BlattWeisskopf);

	virtual std::complex<double> evaluateAmp(dataPoint& point);

	virtual std::shared_ptr<FunctionTree> setupTree(
			allMasses& theMasses,allMasses& toyPhspSample,std::string suffix, ParameterList& params);
protected:
	std::shared_ptr<DoubleParameter> _width;
	double tmp_width;
};

class BreitWignerConf : public basicConf
{
public:
	BreitWignerConf(const boost::property_tree::ptree &pt_) : basicConf(pt_){
		m_mass= pt_.get<double>("mass");
		m_mass_fix= pt_.get<bool>("mass_fix");
		m_mass_min= pt_.get<double>("mass_min");
		m_mass_max= pt_.get<double>("mass_max");
		m_width= pt_.get<double>("width");
		m_width_fix= pt_.get<bool>("width_fix");
		m_width_min= pt_.get<double>("width_min");
		m_width_max= pt_.get<double>("width_max");
		m_mesonRadius= pt_.get<double>("mesonRadius");
		m_spin= pt_.get<unsigned int>("spin");
		m_m= pt_.get<unsigned int>("m");
		m_n= pt_.get<unsigned int>("n");
		m_daughterA= pt_.get<unsigned int>("daughterA");
		m_daughterB= pt_.get<unsigned int>("daughterB");
	}
	virtual void put(boost::property_tree::ptree &pt_){
		basicConf::put(pt_);
		pt_.put("mass", m_mass);
		pt_.put("mass_fix", m_mass_fix);
		pt_.put("mass_min", m_mass_min);
		pt_.put("mass_max", m_mass_max);
		pt_.put("width", m_width);
		pt_.put("width_fix", m_width_fix);
		pt_.put("width_min", m_width_min);
		pt_.put("width_max", m_width_max);
		pt_.put("mesonRadius", m_mesonRadius);
		pt_.put("spin", m_spin);
		pt_.put("m", m_m);
		pt_.put("n", m_n);
		pt_.put("daughterA", m_daughterA);
		pt_.put("daughterB", m_daughterB);
	}
	virtual void update(ParameterList par){
		if(!m_enable) return;
		basicConf::update(par);
		try{// only update parameters if they are found in list
			m_mass= par.GetDoubleParameter("m0_"+m_name)->GetValue();
			m_width= par.GetDoubleParameter("width_"+m_name)->GetValue();
		} catch (BadParameter b) { }
	}

	double m_mass;
	bool m_mass_fix;
	double m_mass_min;
	double m_mass_max;
	double m_width;
	bool m_width_fix;
	double m_width_min;
	double m_width_max;

	double m_mesonRadius;
	int m_spin;
	int m_m;
	int m_n;

	unsigned int m_daughterA;
	unsigned int m_daughterB;
};

class BreitWignerStrategy : public Strategy {
public:
	BreitWignerStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}
	virtual const std::string to_str() const {
		return ("relativistic BreitWigner of "+name);
	}
	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) ;

protected:
	std::string name;
};

class BreitWignerPhspStrategy : public BreitWignerStrategy {
public:
	BreitWignerPhspStrategy(const std::string resonanceName, ParType in):BreitWignerStrategy(resonanceName,in){}
	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out);
};


#endif
