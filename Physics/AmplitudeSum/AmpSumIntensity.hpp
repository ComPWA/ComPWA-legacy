//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     	Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding flatte type resonance, removing root dependence
//-------------------------------------------------------------------------------

#ifndef _AMPSUMINTENSITY_HPP
#define _AMPSUMINTENSITY_HPP

#include <vector>
#include <memory>
#include <map>
#include <string>

#include <boost/property_tree/xml_parser.hpp>

#include "Core/Amplitude.hpp"
#include "Core/Resonance.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Efficiency.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Generator.hpp"

#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

class AmpSumIntensity : public Amplitude {

public:
	AmpSumIntensity( std::string name="", normStyle ns=normStyle::none,
			std::shared_ptr<Efficiency> eff=
					std::shared_ptr<Efficiency>(new UnitEfficiency),
					unsigned int nCalls=30000);

	//! Copy constructor
	AmpSumIntensity( const AmpSumIntensity& copy );

	//! Destructor
	virtual ~AmpSumIntensity(){ /* nothing */ };

	//! Clone function
	virtual AmpSumIntensity* Clone(std::string newName="") const;

	//===================== OPERATORS ==========================
	/** Operator for coherent addition of amplitudes
	 *
	 * @param other
	 * @return
	 */
	const AmpSumIntensity operator+(const AmpSumIntensity& other) const;

	/** Operator for coherent addition of amplitudes
	 *
	 * @param rhs
	 * @return
	 */
	AmpSumIntensity& operator+=(const AmpSumIntensity& rhs) ;

	//===================== LOAD/SAVE CONFIG ====================
	//! Configure resonance from ptree
	virtual void Configure(const boost::property_tree::ptree &pt);

	//! Save resonance from to ptree
	virtual void Save(std::string fileName);

	//!Set prefactor
	virtual void SetPrefactor(std::complex<double> pre);

	//! Get Monte-Carlo precision (number of calls to integration algorithms)
	unsigned int GetMcPrecision() { return _nCalls;}

	//! get maximum value of amplitude with current parameters
	virtual double GetMaxVal( std::shared_ptr<Generator> gen);

	//! Average width of all resonances
	virtual double averageWidth();

	//==================== PRINTING ============================
	//! Print overview over all amplitudes
	virtual void to_str();
	//! print all fit fractions; fitting errors are not available here
	virtual void printFractions();

	//==================== NORMALIZATION/INTEGRATION ==============
	//! Get normalization factor
	virtual const double GetNormalization();

	//! Get amplitdude integral. Result is not efficiency corrected.
	virtual const double GetIntegral();

	//! Get amplitude integral for a set of (sub-) resonances
	virtual const double GetIntegral(std::vector<resonanceItr> resoList);

	//! Calculate inteference integral of two subresonances of @A and @B
	virtual const double GetIntegralInterference(resonanceItr A, resonanceItr B);

	//! calculate integral for a list of resonances
	static const double GetIntegralInterference(
			std::vector<resonanceItr> resList, unsigned int nCalls);

	/** Calculate partial integral over amplitude
	 *
	 * Currently only integration over m23sq and m13sq is supported
	 * @param var1 first integration variables, choose m23sq or m13sq
	 * @param min1 min of first integration variable
	 * @param max1 min of first integration variable
	 * @param var2 second integration variables, choose m23sq or m13sq
	 * @param min2 min of second integration variable
	 * @param max2 max of second integration variable
	 * @return
	 */
	virtual double GetIntValue(std::string var1, double min1, double max1,
			std::string var2, double min2=0, double max2=0);

	/** Calculation of integral
	 *
	 * @param resList List of resonances to integrate
	 * @param eff Efficiency correction (disable -> NULLPTR)
	 * @param nCalls Monte-Carlo precision (number of calls)
	 * @return
	 */
	static double integral(std::vector<resonanceItr> resList,
			std::shared_ptr<Efficiency> eff, int nCalls=30000);

	//==================== EVALUATION =======================
	/**! Evaluate total amplitude
	 * Using current set of parameters at phsp point \point
	 * @param point Point in phase-space
	 * @return Complex function value
	 */
	virtual const std::complex<double>& evaluate( dataPoint& point );

	/**! Evaluate total amplitude
	 * Using current set of parameters at phsp point \point. Amplitude is
	 * multiplied with efficiency of datapoint.
	 */
	virtual const ParameterList& intensity( dataPoint& point );

	/**! Evaluate total amplitude
	 * Using current set of parameters at phsp point \point. Amplitude is
	 * multiplied with efficiency of datapoint.
	 */
	virtual const ParameterList& intensity( std::vector<double> point );

	/**! Evaluate total amplitude
	 * Using current set of parameters at phsp point \point.
	 * No efficiency correction.
	 */
	virtual const ParameterList& intensityNoEff( dataPoint& point );

	virtual const double sliceIntensity(dataPoint& dataP,
			std::complex<double>* reso, unsigned int nResos);

	//==================== FIT FRACTIONS =======================
	virtual void GetFitFractions(ParameterList& parList);

	static void GetFitFractions(ParameterList& parList,
			const Amplitude* amp);

	//================== ACCESS to resonances ====================
	//! Get ID of resonance from name
	virtual int GetIdOfResonance(std::string name);

	//! Get resonance name from ID
	virtual std::string GetNameOfResonance(unsigned int id);

	//! get resonance by @param name
	virtual std::shared_ptr<Resonance> GetResonance(std::string name);

	//! get resonance by @param id
	virtual std::shared_ptr<Resonance> GetResonance(unsigned int id);

	//! List of resonances (enabled AND disabled)
	virtual std::vector<std::shared_ptr<Resonance> >& GetFullListOfResonances(){
		return resoList;
	}

	//! Iterator on first resonance (which is enabled)
	virtual resonanceItr GetResonanceItrFirst(){
		return resonanceItr(
				_resEnabled, resoList.begin(), resoList.end()
		);
	}

	//! Iterator on last resonance (which is enabled)
	virtual resonanceItr GetResonanceItrLast(){
		return resonanceItr(
				_resEnabled, resoList.end(), resoList.end()
		);
	}

	//! List of resonance iterators
	virtual const std::vector<resonanceItr> GetResonanceItrList(){
		std::vector<resonanceItr> resItrList;
		auto itr = GetResonanceItrFirst();
		for( ; itr != GetResonanceItrLast(); ++itr)
			resItrList.push_back(itr);
		return resItrList;
	}

	//========== FunctionTree =============
	//! Check of tree is available
	virtual bool hasTree(){	return 1; }
	//! Getter function for function tree
	virtual std::shared_ptr<FunctionTree> GetTree(
			ParameterList& sample, ParameterList& phspSample,
			ParameterList& toySample);

protected:
	//! Maximum value of amplitude. Necessary for event generation.
	double _maxFcnVal;
	//! Is amplitude maximum already calculated?
	bool _calcMaxFcnVal;
	//! calculate maximum value of amplitude with current parameters
	virtual void calcMaxVal( std::shared_ptr<Generator> gen);
	//! List of resonances
	std::vector<std::shared_ptr<Resonance> > resoList;
	//! Type of normalization
	normStyle _normStyle;
	//! precision for numeric integration
	unsigned int _nCalls;

	//========== FunctionTree =============
	/**Setup Basic Tree
	 *
	 * @param theMasses data sample
	 * @param toyPhspSample sample of flat toy MC events for normalization of the resonances
	 * @param opt Which tree should be created? "data" data Tree, "norm" normalization tree
	 * with efficiency corrected toy phsp sample or "normAcc" normalization tree with sample
	 * of accepted flat phsp events
	 */
	std::shared_ptr<FunctionTree> setupBasicTree(ParameterList& sample,
			ParameterList& phspSample);
};

#endif
