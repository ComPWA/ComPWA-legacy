//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding correct g1s
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef AMP_FLATTE_RES
#define AMP_FLATTE_RES

#include <vector>

#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpKinematics.hpp"
#include "Physics/AmplitudeSum/AmpWigner.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Physics/AmplitudeSum/NonResonant.hpp"

using namespace std;

class AmpFlatteRes : public AmpAbsDynamicalFunction, public AmpKinematics {
public:
	AmpFlatteRes(const char *name,
			std::shared_ptr<DoubleParameter> resMass,
			std::shared_ptr<DoubleParameter> mesonRadius,
			std::shared_ptr<DoubleParameter> motherRadius,
			std::shared_ptr<DoubleParameter> g1, std::shared_ptr<DoubleParameter> g2,
			double _g2_partA, double _g2_partB,
			int _subsys, int resSpin, int m, int n, double resRadius) ;
	virtual ~AmpFlatteRes();

	/** Dynamical function for two coupled channel approach
	 *
	 * @param mSq center-of-mass energy^2 (=s)
	 * @param mR mass of resonaces
	 * @param massA1 mass of first particle of signal channel
	 * @param massA2 mass of second particle of signal channel
	 * @param gA coupling constant for signal channel
	 * @param massB1 mass of first particle of second channel
	 * @param massB2 mass of second particle of second channel
	 * @param gB coupling constant for secondchannel
	 * @param J resonace spin
	 * @param mesonRadius 1/interaction length (needed for barrier factors)
	 * @return
	 */
	static std::complex<double> dynamicalFunction(double mSq, double mR,
			double massA1, double massA2, double gA,
			double massB1, double massB2, double gB,
			unsigned int J, double mesonRadius);

	virtual void initialise() { };
	std::complex<double> evaluate(dataPoint& point) { return ( _norm*evaluateAmp(point)*evaluateWignerD(point) ); }
	virtual std::complex<double> evaluateAmp(dataPoint& point) ;

	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	unsigned int getNParams(){ return nParams;}

protected:
	unsigned int nParams;
	double _g2_partA;//hidden channel: mass particle A
	double _g2_partB; //hidden channel: mass particle B
	std::shared_ptr<DoubleParameter> _g2, _g1;
	bool foundMasses;
	unsigned int id23, id13;
	double mesonRadius;
};

class FlatteConf : public basicConf
{
public:
	virtual ~FlatteConf() { }
	FlatteConf(const boost::property_tree::ptree &pt_) : basicConf(pt_){
		m_mass= pt_.get<double>("mass");
		m_mass_fix= pt_.get<bool>("mass_fix");
		m_mass_min= pt_.get<double>("mass_min");
		m_mass_max= pt_.get<double>("mass_max");
		m_mesonRadius= pt_.get<double>("mesonRadius");
		m_spin= pt_.get<unsigned int>("spin");
		m_m= pt_.get<unsigned int>("m");
		m_n= pt_.get<unsigned int>("n");
		m_daughterA= pt_.get<unsigned int>("daughterA");
		m_daughterB= pt_.get<unsigned int>("daughterB");
		m_g1= pt_.get<double>("g1");
		m_g1_fix= pt_.get<bool>("g1_fix");
		m_g1_min= pt_.get<double>("g1_min");
		m_g1_max= pt_.get<double>("g1_max");
		m_g2= pt_.get<double>("g2");
		m_g2_part1= pt_.get<std::string>("g2_part1");
		m_g2_part2= pt_.get<std::string>("g2_part2");

	}
	virtual void put(boost::property_tree::ptree &pt_){
		basicConf::put(pt_);
		pt_.put("mass", m_mass);
		pt_.put("mass_fix", m_mass_fix);
		pt_.put("mass_min", m_mass_min);
		pt_.put("mass_max", m_mass_max);
		pt_.put("mesonRadius", m_mesonRadius);
		pt_.put("spin", m_spin);
		pt_.put("m", m_m);
		pt_.put("n", m_n);
		pt_.put("daughterA", m_daughterA);
		pt_.put("daughterB", m_daughterB);
		pt_.put("g1", m_g1);
		pt_.put("g1_fix", m_g1_fix);
		pt_.put("g1_min", m_g1_min);
		pt_.put("g1_max", m_g1_max);
		pt_.put("g2", m_g2);
		pt_.put("g2_part1", m_g2_part1);
		pt_.put("g2_part2", m_g2_part2);
	}
	virtual void update(ParameterList par){
		basicConf::update(par);
		try{// only update parameters if they are found in list
		m_mass= par.GetDoubleParameter("m0_"+m_name)->GetValue();
		m_g1= par.GetDoubleParameter("g1_"+m_name)->GetValue();
		m_g2= par.GetDoubleParameter("g2_"+m_name)->GetValue();
		} catch (BadParameter b) { }
	}
	double m_mass;
	bool m_mass_fix;
	double m_mass_min;
	double m_mass_max;

	double m_mesonRadius;
	unsigned int m_spin;
	unsigned int m_m;
	unsigned int m_n;

	unsigned int m_daughterA; //TODO: better reference
	unsigned int m_daughterB; //TODO: better reference
	double m_g1;
	double m_g1_fix;
	double m_g1_min;
	double m_g1_max;
	double m_g2;
	std::string m_g2_part1;
	std::string m_g2_part2;
};

class FlatteStrategy : public Strategy {
public:
	FlatteStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}
	virtual const std::string to_str() const { return ("flatte amplitude of "+name); }
	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out);

protected:
	std::string name;
};

class FlattePhspStrategy : public Strategy {
public:
	FlattePhspStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}
	virtual const std::string to_str() const { return ("flatte amplitude of "+name); }
	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out);

protected:
	std::string name;
};
#endif
