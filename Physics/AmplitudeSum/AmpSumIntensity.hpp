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

#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

class AmpSumIntensity : public Amplitude {

public:

	AmpSumIntensity(AmplitudeSetup ini, normStyle ns,
			std::shared_ptr<Efficiency> eff, unsigned int nCalls);
	AmpSumIntensity(AmplitudeSetup ini,
			std::shared_ptr<Efficiency> eff, unsigned int nCalls);
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

	//! setting new parameterList
	virtual void setParameterList(ParameterList& par);
	//! evaluate total amplitude using parameters \par at phsp point \point
	virtual const ParameterList& intensity(dataPoint& point, ParameterList& par);
	//! evaluate total amplitude using current set of parameters at phsp point \point. Amplitude is multiplied with efficiency of datapoint.
	virtual const ParameterList& intensity(dataPoint& point);
	//! evaluate total amplitude using current set of parameters at phsp point \point. No efficiency correction.
	virtual const ParameterList& intensityNoEff(dataPoint& point);
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
	virtual double getIntValue(std::string var1, double min1, double max1,
			std::string var2, double min2=0, double max2=0);

	//---------- get resonance parameters -------------
	//! \return Number of resonances
	virtual unsigned int GetNumberOfResonances();
	//! convert resonance \param name to id
	int GetIdOfResonance(std::string name);
	//! convert resonance \param id to name
	std::string GetNameOfResonance(unsigned int id);
	//! get magnitude of resonance \param name
	virtual double GetMagnitude(std::string name);
	//! get magnitude of resonance \param ID
	virtual double GetMagnitude(unsigned int id);
	//! get phase of resonance \param name
	virtual double GetPhase(std::string name);
	//! get phase of resonance \param ID
	virtual double GetPhase(unsigned int id);
	//! get total integral for resonance \param name
	virtual double GetIntegral(std::string name);
	//! get total integral for resonance \param id
	virtual double GetIntegral(unsigned int id);
	//! get fit fraction for resonance \param name
	virtual double GetFraction(std::string name);
	//! get fit fraction for resonance \param id
	virtual double GetFraction(unsigned int id);
	//! get resonance by @param name
	virtual std::shared_ptr<AmpAbsDynamicalFunction> GetResonance(std::string name);
	//! get resonance by @param id
	virtual std::shared_ptr<AmpAbsDynamicalFunction> GetResonance(unsigned int id);
	/*!Get AmplitudeSetup
	 * AmpltidueSetup object is updated with current parameters and a pointer is returned.
	 *
	 * @return AmplitudeSetup
	 */
	AmplitudeSetup* GetAmplitudeSetup();

	//---------- related to FunctionTree -------------
	//! Check of tree is available
	virtual bool hasTree(){	return 1; }
	//! Getter function for function tree
	virtual std::shared_ptr<FunctionTree> GetTree(allMasses& theMasses,allMasses& toyPhspSample,
			std::string suffix=""){
		return setupBasicTree(theMasses,toyPhspSample, suffix);
	}

protected:
	//! Initialize amplitude sum from amplitude setup
	void init();
	//! Maximum value of amplitude. Necessary for event generation.
	double _maxFcnVal;
	//! Is amplitude maximum already calculated?
	bool _calcMaxFcnVal;
	//! calculate maximum value of amplitude with parameters \par
	virtual void calcMaxVal(ParameterList& par ,std::shared_ptr<Generator> gen);
	//! calculate maximum value of amplitude with current parameters
	virtual void calcMaxVal( std::shared_ptr<Generator> gen);
	//! Efficiency object
	std::shared_ptr<Efficiency> eff_;
	//! List of resonances
	std::vector<std::shared_ptr<AmpAbsDynamicalFunction> > resoList;
	//! Amplitude setup
	AmplitudeSetup ampSetup;
	//! Type of normalization
	normStyle _normStyle;
	//! List of parameters
	ParameterList params;
	//! Update amplitude setup with current parameters
	void updateAmplitudeSetup();
	//! precision for numeric integration
	unsigned int _nCalls;

	//---------- related to FunctionTree -------------
	/**Setup Basic Tree
	 *
	 * @param theMasses data sample
	 * @param toyPhspSample sample of flat toy MC events for normalization of the resonances
	 * @param opt Which tree should be created? "data" data Tree, "norm" normalization tree
	 * with efficiency corrected toy phsp sample or "normAcc" normalization tree with sample
	 * of accepted flat phsp events
	 */
	std::shared_ptr<FunctionTree> setupBasicTree(allMasses& theMasses,allMasses& toyPhspSample,
			std::string suffix="");

private:


};

#endif
