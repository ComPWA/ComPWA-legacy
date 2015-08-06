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

#ifndef AMP_FLATTE3CH_RES
#define AMP_FLATTE3CH_RES

#include <vector>
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"

using namespace std;

class AmpFlatteRes3Ch : public AmpAbsDynamicalFunction, public AmpKinematics {
public:
	AmpFlatteRes3Ch(const char *name,
			std::shared_ptr<DoubleParameter> resMass,
			std::shared_ptr<DoubleParameter> mesonRadius,
			std::shared_ptr<DoubleParameter> motherRadius,
			std::shared_ptr<DoubleParameter> g1,
			std::shared_ptr<DoubleParameter> g2, double _g2_partA, double _g2_partB,
			std::shared_ptr<DoubleParameter> g3, double _g3_partA, double _g3_partB,
			int _subsys, int resSpin, int m, int n, double resRadius) ;
	virtual ~AmpFlatteRes3Ch();

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
			double massC1, double massC2, double gC,
			unsigned int J, double mesonRadius);

	virtual void initialise() { };
	std::complex<double> evaluate(const dataPoint& point) { return ( _norm*evaluateAmp(point)*evaluateWignerD(point) ); }
	virtual std::complex<double> evaluateAmp(const dataPoint& point) ;

	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	unsigned int getNParams(){ return nParams;}

protected:
	unsigned int nParams;
	double _g2_partA;//hidden channel: mass particle A
	double _g2_partB; //hidden channel: mass particle B
	double _g3_partA;//hidden channel: mass particle A
	double _g3_partB; //hidden channel: mass particle B
	std::shared_ptr<DoubleParameter> _g3, _g2, _g1;
	bool foundMasses;
	unsigned int id23, id13;
	double mesonRadius;
};

class Flatte3ChConf : public FlatteConf
{
public:
	virtual ~Flatte3ChConf() { }
	Flatte3ChConf(const boost::property_tree::ptree &pt_);
	virtual void put(boost::property_tree::ptree &pt_);
	virtual void update(ParameterList par);

	double m_g3;
	std::string m_g3_part1;
	std::string m_g3_part2;
};


class Flatte3ChStrategy : public Strategy {
public:
	Flatte3ChStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}
	virtual const std::string to_str() const { return ("flatte3Ch amplitude of "+name); }
	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out);

protected:
	std::string name;
};

class Flatte3ChPhspStrategy : public Strategy {
public:
	Flatte3ChPhspStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}
	virtual const std::string to_str() const { return ("flatte3Ch amplitude of "+name); }
	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out);

protected:
	std::string name;
};
#endif
