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
#include <cmath>
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"

using namespace std;

class AmpFlatteRes : public AmpAbsDynamicalFunction
{
public:

	AmpFlatteRes( normStyle nS=normStyle::one, int calls=30000 );

	AmpFlatteRes(const char *name,
			unsigned int varIdA, unsigned int varIdB,
			std::shared_ptr<DoubleParameter> mag,
			std::shared_ptr<DoubleParameter> phase,
			std::shared_ptr<DoubleParameter> mass,
			Spin spin, Spin m, Spin n, int P, int C,
			std::string mother, std::string particleA, std::string particleB,
			std::shared_ptr<DoubleParameter> mesonRadius,
			std::shared_ptr<DoubleParameter> motherRadius,
			std::shared_ptr<DoubleParameter> g1,
			std::shared_ptr<DoubleParameter> g2, std::string g2_idA, std::string g2_idB,
			std::shared_ptr<DoubleParameter> g3, std::string g3_idA, std::string g3_idB,
			formFactorType type = formFactorType::CrystalBarrel,
			int nCalls=30000, normStyle nS=normStyle::one) ;

	virtual ~AmpFlatteRes();

	//! Clone function
	virtual AmpFlatteRes* Clone(std::string newName="") const{
		auto tmp = (new AmpFlatteRes(*this));
		if(newName != "")
			tmp->SetName(newName);
		return tmp;
	}

	//! Configure resonance from ptree
	virtual void Configure(boost::property_tree::ptree::value_type const& v,
			ParameterList& list);

	//! Save resonance from to ptree
	virtual void Save(boost::property_tree::ptree&);

	//! Check of parameters have changed and normalization has to be recalculatecd
	virtual void CheckModified();

	//! Print resonance parameters
	std::string to_str() const;

	//! Calculation integral |dynamical amplitude|^2
	virtual double GetIntegral();

	//! Get resonance width
	virtual double GetWidth() const {
		return std::abs( couplingToWidth(_mass->GetValue(),_mass->GetValue(), _g1->GetValue(),
			_mass1, _mass2, _spin, _mesonRadius->GetValue(), _ffType) );
	}

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
			double massC1, double massC2, double gC,
			unsigned int J, double mesonRadius, formFactorType ffType=formFactorType::CrystalBarrel );

	virtual std::complex<double> EvaluateAmp(dataPoint& point) ;

	virtual std::shared_ptr<FunctionTree> SetupTree(
			ParameterList& sample, ParameterList& toySample,std::string suffix);

protected:
	//Initialize masses
	void initialize();

	double _g2_massA, _g2_massB, _g3_massA, _g3_massB;
	std::string _g2_idA, _g2_idB, _g3_idA, _g3_idB;
	std::shared_ptr<DoubleParameter> _g3, _g2, _g1;
	double tmp_g3, tmp_g2, tmp_g1;
	double _g1_writeByName, _g2_writeByName, _g3_writeByName;
};

class FlatteStrategy : public Strategy
{
public:
	FlatteStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}
	virtual const std::string to_str() const { return ("flatte amplitude of "+name); }
	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out);

protected:
	std::string name;
};

#endif
