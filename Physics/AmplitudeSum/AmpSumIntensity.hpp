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

#include "Core/Amplitude.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Efficiency.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Generator.hpp"

#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

class AmpSumIntensity : public Amplitude {

public:
	enum normStyle {
		none, /*!< no normaliztion between Amplitudes. */
		/*!< all amplitudes are normalized to one.
		 *  The normalization factor is \f$ 1/\sqrt(\int |A|^2)\f$ */
		one

	};
	AmpSumIntensity(AmplitudeSetup ini, normStyle ns,
			std::shared_ptr<Efficiency> eff, unsigned int nCalls);
	AmpSumIntensity(AmplitudeSetup ini,
			std::shared_ptr<Efficiency> eff, unsigned int nCalls);
//	AmpSumIntensity(const AmpSumIntensity& other);
	//! Destructor
	virtual ~AmpSumIntensity(){ /* nothing */ };
	//! Clone function
	virtual AmpSumIntensity* Clone(){
		return (new AmpSumIntensity(*this));
	}

	//! set efficiency
	virtual void setEfficiency(std::shared_ptr<Efficiency> eff) { eff_ = eff; };
	//! wrapper function for function value times efficiency at point x
	double evaluateEff(double x[], size_t dim);
	//! wrapper function for function value at point x
	double evaluate(double x[], size_t dim);
	//! normalization integral for parameters \par (doesn't include calculated efficiency)
	virtual const double normalization(ParameterList& par);
	//! normalization integral (doesn't include calculated efficiency)
	virtual const double normalization();
	//! normalization integral for parameters \par (includes calculated efficiency)
	virtual const double integral(ParameterList& par);
	//! normalization integral (includes calculated efficiency)
	virtual const double integral();
	//! get maximum value of amplitude with current parameters
	virtual double getMaxVal( std::shared_ptr<Generator> gen);
	//! get maximum value of amplitude with parameters \par
	virtual double getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen);

	virtual std::complex<double> getFirstAmp(dataPoint& point, ParameterList& par);
	virtual std::complex<double> getFirstReso(dataPoint& point, ParameterList& par);
	virtual std::complex<double> getFirstBW(dataPoint& point, ParameterList& par);

	//! setting new parameterList
	virtual void setParameterList(ParameterList& par);
	//! evaluate total amplitude using parameters \par at phsp point \point
	virtual const ParameterList& intensity(const dataPoint& point, ParameterList& par);
	//! evaluate total amplitude using current set of parameters at phsp point \point. Amplitude is multiplied with efficiency of datapoint.
	virtual const ParameterList& intensity(const dataPoint& point);
	//! evaluate total amplitude using current set of parameters at phsp point \point. No efficiency correction.
	virtual const ParameterList& intensityNoEff(const dataPoint& point);
	//! evaluate total amplitude using current set of parameters at phsp point \point. Amplitude is multiplied with efficiency of datapoint.
	virtual const ParameterList& intensity(std::vector<double> point, ParameterList& par);
	virtual const double sliceIntensity(dataPoint& dataP, ParameterList& par,std::complex<double>* reso, unsigned int nResos);

	//! fill internal parameter list with (start) parameter
	virtual bool copyParameterList(ParameterList& outPar);
	//! print overview over all amplitudes
	virtual void printAmps();
	//! print all fit fractions; fitting errors are not available here
	virtual void printFractions();
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
	virtual double getIntValue(std::string var1, double min1, double max1, std::string var2, double min2=0, double max2=0);

	//---------- get resonance parameters -------------
	//! \return Number of resonances
	virtual unsigned int getNumberOfResonances() { return totAmp.getNumberOfResonances(); }
	//! convert resonance \param name to id
	int getIdOfResonance(std::string name){ return totAmp.getAmpId(name); }
	//! convert resonance \param id to name
	std::string getNameOfResonance(unsigned int id){ return totAmp.getAmpName(id); }
	//! get magnitude of resonance \param name
	virtual double getAmpMagnitude(std::string name) { return totAmp.getAmpMagnitude(name); };
	//! get magnitude of resonance \param ID
	virtual double getAmpMagnitude(unsigned int id) { return totAmp.getAmpMagnitude(id); };
	//! get phase of resonance \param name
	virtual double getAmpPhase(std::string name) { return totAmp.getAmpPhase(name); };
	//! get phase of resonance \param ID
	virtual double getAmpPhase(unsigned int id) { return totAmp.getAmpPhase(id); };
	//! get total integral for resonance \param name
	virtual double getAmpIntegral(std::string name) { return totAmp.getAmpIntegral(name); };
	//! get total integral for resonance \param id
	virtual double getAmpIntegral(unsigned int id) { return totAmp.getAmpIntegral(id); };
	//! get fit fraction for resonance \param name
	virtual double getAmpFraction(std::string name) {
		return totAmp.getAmpStrength(name)/integral();
	};
	//! get fit fraction for resonance \param id
	virtual double getAmpFraction(unsigned int id) {
		return totAmp.getAmpStrength(id)/integral();
	};
	//! get resonance by @param name
	virtual std::shared_ptr<AmpAbsDynamicalFunction> getResonance(std::string name) {
		return totAmp.getResonance(name);
	};
	//! get resonance by @param id
	virtual std::shared_ptr<AmpAbsDynamicalFunction> getResonance(unsigned int id) {
		return totAmp.getResonance(id);
	};
	//! Get pointer to amplitude sum
	AmpSumOfAmplitudes* getAmpSum(){ return &totAmp; }


	/*!Get AmplitudeSetup
	 * AmpltidueSetup object is updated with current parameters and a pointer is returned.
	 *
	 * @return AmplitudeSetup
	 */
	AmplitudeSetup* GetAmplitudeSetup() {
		updateAmplitudeSetup();
		return &ampSetup;
	}

	//---------- related to FunctionTree -------------
	//! Check of tree is available
	virtual bool hasTree(){	return 1; }
	//! Getter function for function tree
	virtual std::shared_ptr<FunctionTree> getAmpTree(allMasses& theMasses,allMasses& toyPhspSample, std::string suffix=""){
		return setupBasicTree(theMasses,toyPhspSample, suffix);
	}

protected:
	void init();
	//! calculate maximum value of amplitude with parameters \par
	virtual void calcMaxVal(ParameterList& par ,std::shared_ptr<Generator> gen);
	//! calculate maximum value of amplitude with current parameters
	virtual void calcMaxVal( std::shared_ptr<Generator> gen);
	//! calculate normalization of resonance \param amp
	double normReso(std::shared_ptr<AmpAbsDynamicalFunction> amp);

	std::shared_ptr<Efficiency> eff_;
	bool _calcMaxFcnVal;
	bool _calcNorm;
	double _maxFcnVal;
	AmpSumOfAmplitudes totAmp;
	AmplitudeSetup ampSetup;

	double maxVal;
	normStyle _normStyle;
	ParameterList params;

	void updateAmplitudeSetup();
	unsigned int _nCalls; //! precision for numeric integration

	//---------- related to FunctionTree -------------
	/**Setup Basic Tree
	 *
	 * @param theMasses data sample
	 * @param toyPhspSample sample of flat toy MC events for normalization of the resonances
	 * @param opt Which tree should be created? "data" data Tree, "norm" normalization tree
	 * with efficiency corrected toy phsp sample or "normAcc" normalization tree with sample
	 * of accepted flat phsp events
	 */
	std::shared_ptr<FunctionTree> setupBasicTree(allMasses& theMasses,allMasses& toyPhspSample, std::string suffix="");

private:

};

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
