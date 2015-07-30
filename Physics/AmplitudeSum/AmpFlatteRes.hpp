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
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Physics/AmplitudeSum/NonResonant.hpp"

using namespace std;


class AmpFlatteRes : public AmpAbsDynamicalFunction {
public:
	AmpFlatteRes(const char *name,
			std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
			std::shared_ptr<DoubleParameter> mass, int subSys, Spin spin, Spin m, Spin n,
			std::shared_ptr<DoubleParameter> mesonRadius,
			std::shared_ptr<DoubleParameter> motherRadius,
			std::shared_ptr<DoubleParameter> g1, std::shared_ptr<DoubleParameter> g2,
			double g2_partA, double g2_partB, int nCalls=30000, normStyle nS=normStyle::one) ;
	virtual ~AmpFlatteRes();

	/** Dynamical function for two coupled channel approach
	 *
	 * @param mSq center-of-mass energy^2 (=s)
	 * @param mR mass of resonances
	 * @param massA1 mass of first particle of signal channel
	 * @param massA2 mass of second particle of signal channel
	 * @param gA coupling constant for signal channel
	 * @param massB1 mass of first particle of second channel
	 * @param massB2 mass of second particle of second channel
	 * @param gB coupling constant for second channel
	 * @param J resonance spin
	 * @param mesonRadius 1/interaction length (needed for barrier factors)
	 * @return
	 */
	static std::complex<double> dynamicalFunction(double mSq, double mR,
			double massA1, double massA2, double gA,
			double massB1, double massB2, double gB,
			unsigned int J, double mesonRadius);

	virtual std::complex<double> evaluateAmp(dataPoint& point) ;

	virtual std::shared_ptr<FunctionTree> setupTree(
			allMasses& theMasses,allMasses& toyPhspSample,std::string suffix, ParameterList& params);
protected:
	double _g2_partA;//hidden channel: mass particle A
	double _g2_partB; //hidden channel: mass particle B
	std::shared_ptr<DoubleParameter> _g2, _g1;
	double tmp_g2, tmp_g1, tmp_mass;
	double mesonRadius;
};

class FlatteConf : public basicConf
{
public:
	virtual ~FlatteConf() { }
	FlatteConf(const boost::property_tree::ptree &pt_);
	virtual void put(boost::property_tree::ptree &pt_);
	virtual void update(ParameterList par);

	double m_mass;
	bool m_mass_fix;
	double m_mass_min;
	double m_mass_max;

	double m_mesonRadius;
	int m_spin;
	int m_m;
	int m_n;

	unsigned int m_daughterA; //TODO: better reference
	unsigned int m_daughterB; //TODO: better reference
	double m_g1;
	double m_g1_fix;
	double m_g1_min;
	double m_g1_max;
	double m_g2;
	double m_g2_fix;
	double m_g2_min;
	double m_g2_max;
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
