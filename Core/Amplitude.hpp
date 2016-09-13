//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding UnitAmp
//-------------------------------------------------------------------------------
//! Physics Interface Base-Class.
/*! \class Amplitude
 * @file Amplitude.hpp
 * This class provides the interface to the model which tries to describe the
 * intensity. As it is pure virtual, one needs at least one implementation to
 * provide an model for the analysis which calculates intensities for an event on
 * basis model parameters. If a new physics-model is derived from and fulfills
 * this base-class, no change in other modules are necessary to work with the new
 * physics module.
 */

#ifndef PIFBASE_HPP_
#define PIFBASE_HPP_

#include <vector>
#include <memory>
#include <math.h>

#include "Core/Resonance.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "DataReader/Data.hpp"

#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Generator.hpp"


class Amplitude
{

public:

	Amplitude(std::string name="",
			std::shared_ptr<Efficiency> eff	= std::shared_ptr<Efficiency>(new UnitEfficiency)) :
				_name(name), eff_(eff) { }

	virtual ~Amplitude() { /* nothing */	}

	virtual Amplitude* Clone(std::string newName="") const = 0;

	//============ SET/GET =================
	//! Get name of amplitude
	virtual std::string GetName() const { return _name; }

	//! Set name of amplitude
	virtual void SetName(std::string name) { _name = name; }

	//! Get efficiency
	virtual std::shared_ptr<Efficiency> GetEfficiency() { return eff_; };

	//! Set efficiency
	virtual void SetEfficiency(std::shared_ptr<Efficiency> eff) { eff_ = eff; };

	/**! Set efficiencies for a vector of amplitudes
	 *
	 * @param ampVec Vector of amplitudes
	 * @param eff New efficiency object
	 */
	static void SetAmpEfficiency( std::vector<std::shared_ptr<Amplitude> > ampVec,
		std::shared_ptr<Efficiency> eff);

	/** Get maximum value of amplitude
	 * Maximum is numerically calculated using a random number generator
	 * @param gen Random number generator
	 * @return
	 */
	virtual double GetMaxVal(std::shared_ptr<Generator> gen) = 0;

	//============= PRINTING =====================
	//! Print amplitude to logging system
	virtual void to_str() = 0;

	//=========== INTEGRATION/NORMALIZATION =================
	/** Calculate normalization of amplitude.
	 * The integral includes efficiency correction
	 */
	virtual const double GetNormalization() = 0;

	/** Calculate integral of amplitude.
	 * The integral does not include efficiency correction
	 */
	virtual const double GetIntegral() = 0;

	/** Calculate integral of amplitudes for a given vector.
	 * The integral does not include efficiency correction
	 *
	 * @param resoList Vector of amplitudes
	 * @return
	 */
	virtual const double GetIntegral(std::vector<resonanceItr> resoList) {
		return 0;
	};

	/** Calculate interference integral
	 *
	 * @param A First resonance
	 * @param B Second resonance
	 * @return a*conj(b)+conj(a)*b
	 */
	virtual const double GetIntegralInterference(resonanceItr A, resonanceItr B)
	{ return -999; };

	/** Integral value of amplitude in certain boundary
	 * Used for plotting a projection of a function in \p var1 in
	 * bin [\p min1, \p min2]. In this case we have to integrate over an
	 * other variable \p var2
	 * @param var1 first variable
	 * @param min1 minimal value of first variable
	 * @param max1 maximal value of first variable
	 * @param var2 second variable
	 * @param min2 minimal value of second variable
	 * @param max2 maximal value of second variable
	 * @return
	 */
	virtual double GetIntValue(std::string var1, double min1, double max1,
			std::string var2, double min2, double max2) = 0;

	//=========== EVALUATION =================
	/** Calculate intensity of amplitude at point in phase-space
	 *
	 * @param point Data point
	 * @return
	 */
	virtual const ParameterList& intensity( dataPoint& point ) = 0;

	/** Calculate intensity of amplitude at point in phase-space
	 *
	 * @param point Data point
	 * @return
	 */
	virtual const ParameterList& intensity(	std::vector<double> point ) = 0;

	/** Calculate intensity of amplitude at point in phase-space
	 * Intensity is calculated excluding efficiency correction
	 * @param point Data point
	 * @return
	 */
	virtual const ParameterList& intensityNoEff( dataPoint& point ) = 0;

	/**! Evaluate interference term of two resonances
	 *
	 * @param point Data point
	 * @param A First resonance
	 * @param B Second resonance
	 * @return
	 */
	static double intensityInterference(dataPoint& point,
			resonanceItr A, resonanceItr B);

	//=========== PARAMETERS =================
	/** Update parameters of resonance
	 *
	 * @param par New list of parameters
	 */
	virtual void UpdateParameters(ParameterList& par);

	/**! Update parameters of vector of amplitudes
	 *
	 * @param ampVec Vector of amplitudes
	 * @param list New list of parameters
	 */
	static void UpdateAmpParameterList(
			std::vector<std::shared_ptr<Amplitude> > ampVec,
			ParameterList& list);

	/**! Add amplitude parameters to list
	 * Add parameters only to list if not already in
	 * @param list Parameter list to be filled
	 */
	virtual void FillParameterList(ParameterList& list) const;

	/**! Add amplitude parameters to list
	 * Add parameter only if not already in list
	 * @param ampVec Vector of amplitudes
	 * @param list Parameter List to be filled
	 */
	static void FillAmpParameterToList(
			std::vector<std::shared_ptr<Amplitude> > ampVec,
			ParameterList& list);

	//! Calculate & fill fit fractions of this amplitude to ParameterList
	virtual void GetFitFractions(ParameterList& parList) = 0;

	/**! Calculate & fill fit fractions of amplitude Vector to list
	 *
	 * @param ampVec Vector of amplitudes
	 * @param list Parameter List to be filled
	 */
	static void GetAmpFitFractions(
			std::vector<std::shared_ptr<Amplitude> > ampVec,
			ParameterList& list);

	//============= ACCESS TO RESONANCES ================
	//! Iterator on first resonance (which is enabled)
	virtual resonanceItr GetResonanceItrFirst() {};

	//! Iterator on last resonance (which is enabled)
	virtual resonanceItr GetResonanceItrLast() {};

	//! Iterator on last resonance (which is enabled)
	virtual const std::vector<resonanceItr> GetResonanceItrList() {};

	//========== FUNCTIONTREE =============
	//! Check of tree is available
	virtual bool hasTree(){ return 0; }

	//! Getter function for basic amp tree
	virtual std::shared_ptr<FunctionTree> GetTree(
			ParameterList&, ParameterList&, ParameterList&) {
		return std::shared_ptr<FunctionTree>();
	}

	/**! Check if amplitudes have a FunctionTree
	 *
	 * @param ampVec
	 * @return
	 */
	static bool AmpHasTree(std::vector<std::shared_ptr<Amplitude> > ampVec);

	unsigned int GetMcPrecision() { return 30000;}

protected:
	//! Name
	std::string _name;

	//! need to store this object for boost::filter_iterator
	resIsEnabled _resEnabled;

	//! Amplitude value
	ParameterList result;

	//! List of interal parameters
	ParameterList params;

	//! Efficiency object
	std::shared_ptr<Efficiency> eff_;
};
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
class GaussAmp : public Amplitude
{
public:
	GaussAmp(const char *name, DoubleParameter _resMass,
			DoubleParameter _resWidth){
		params.AddParameter(
				std::shared_ptr<DoubleParameter>(new DoubleParameter(_resMass))
		);
		params.AddParameter(
				std::shared_ptr<DoubleParameter>(new DoubleParameter(_resWidth))
		);
		initialise();
	}

	GaussAmp(const char *name, double _resMass, double _resWidth){
		params.AddParameter(
				std::shared_ptr<DoubleParameter>(
						new DoubleParameter("mass",_resMass)
				)
		);
		params.AddParameter(
				std::shared_ptr<DoubleParameter>(
						new DoubleParameter("width",_resWidth)
				)
		);
		initialise();
	}

	//! Clone function
	GaussAmp* Clone(std::string newName="") const {
		auto tmp = (new GaussAmp(*this));
		tmp->SetName(newName);
		return tmp;
	}

	virtual void initialise() {
		result.AddParameter(
				std::shared_ptr<DoubleParameter>(
						new DoubleParameter("GaussAmpResult")
				)
		);
		if(Kinematics::instance()->GetNVars()!=1)
			throw std::runtime_error("GaussAmp::initialize() | "
					"this amplitude is for two body decays only!");
	};
	//! Clone function
	virtual GaussAmp* Clone(std::string newName=""){
		auto tmp = (new GaussAmp(*this));
		tmp->SetName(newName);
		return tmp;
	}

	virtual void to_str() { };

	virtual double GetIntValue(std::string var1, double min1, double max1,
			std::string var2, double min2, double max2) { return 0; }

	//! Get integral
	virtual const double GetIntegral(){
		return (params.GetDoubleParameter(1)->GetValue() * std::sqrt(2*M_PI));
	}
	virtual const double GetNormalization() { return GetIntegral(); }

	virtual double GetMaxVal(std::shared_ptr<Generator> gen){
		double mass = params.GetDoubleParameter(0)->GetValue();
		std::vector<double> m; m.push_back(mass*mass);
		dataPoint p(m);
		intensity(p);
		return (result.GetDoubleParameterValue(0));
	}

	virtual const ParameterList& intensity(dataPoint& point) {
		double mass = params.GetDoubleParameter(0)->GetValue();
		double width = params.GetDoubleParameter(1)->GetValue();
		double sqrtS = std::sqrt(point.getVal(0));

		std::complex<double> gaus(
				std::exp(-1*(sqrtS-mass)*(sqrtS-mass)/width/width/2.),
				0
		);

		if(gaus.real() != gaus.real())
			BOOST_LOG_TRIVIAL(error)<<"GaussAmp::intensity() | result NaN!";
		result.SetParameterValue(0,std::norm(gaus));
		return result;
	}

	virtual const ParameterList& intensity(std::vector<double> point){
		dataPoint dataP(point);
		return intensity(dataP);
	}

	virtual const ParameterList& intensityNoEff(dataPoint& point){
		return intensity(point);
	}

	virtual void GetFitFractions(ParameterList& parList) {}

};
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/**! UnitAmp
 *
 * Example implementation of Amplitude with the function value 1.0 at all
 * points in PHSP. It is used to test the likelihood normalization.
 */
class UnitAmp : public Amplitude
{
public:
	UnitAmp() {
		result.AddParameter(
				std::shared_ptr<DoubleParameter>(
						new DoubleParameter("AmpSumResult")
				)
		);
		eff_ = std::shared_ptr<Efficiency>(new UnitEfficiency());
	}

	virtual ~UnitAmp()	{ /* nothing */	}

	virtual UnitAmp* Clone(std::string newName="") const {
		auto tmp = new UnitAmp(*this);
		tmp->SetName(newName);
		return tmp;
	}

	virtual void to_str() { }

	virtual double GetMaxVal(std::shared_ptr<Generator> gen) { return 1; }

	virtual const double GetIntegral() {
		return Kinematics::instance()->GetPhspVolume();
	}
	virtual const double GetNormalization() {
		BOOST_LOG_TRIVIAL(info) << "UnitAmp::normalization() | "
				"normalization not implemented!";
		return 1;
	}
	virtual double GetIntValue(std::string var1, double min1, double max1,
			std::string var2, double min2, double max2) { return 0; }

	virtual const ParameterList& intensity(dataPoint& point) {
		result.SetParameterValue(0,eff_->evaluate(point));
		return result;
	}
	virtual const ParameterList& intensity(std::vector<double> point){
		dataPoint dataP(point);
		return intensity(dataP);
	}
	virtual const ParameterList& intensityNoEff(dataPoint& point) {
		result.SetParameterValue(0,1.0);
		return result;
	}
	virtual void GetFitFractions(ParameterList& parList) {}

	//========== FunctionTree =============
	//! Check of tree is available
	virtual bool hasTree(){ return 1; }

	//! Getter function for basic amp tree
	virtual std::shared_ptr<FunctionTree> getAmpTree(ParameterList& sample,
			ParameterList& toySample, std::string suffix){
		return setupBasicTree(sample,toySample, suffix);
	}

protected:

	/**Setup Basic Tree
	 *
	 * @param theMasses data sample
	 * @param toyPhspSample sample of flat toy MC events for normalization of the resonances
	 * @param opt Which tree should be created? "data" data Tree, "norm" normalization tree
	 * with efficiency corrected toy phsp sample or "normAcc" normalization tree with sample
	 * of accepted flat phsp events
	 */
	std::shared_ptr<FunctionTree> setupBasicTree(ParameterList& sample,
			ParameterList& toySample, std::string suffix) {

		int sampleSize = sample.GetMultiDouble(0)->GetNValues();
		int toySampleSize = toySample.GetMultiDouble(0)->GetNValues();

		BOOST_LOG_TRIVIAL(debug) << "UnitAmp::setupBasicTree() generating new tree!";
		if(sampleSize==0){
			BOOST_LOG_TRIVIAL(error) << "UnitAmp::setupBasicTree() data sample empty!";
			return std::shared_ptr<FunctionTree>();
		}
		std::shared_ptr<FunctionTree> newTree(new FunctionTree());
		//std::shared_ptr<MultAll> mmultDStrat(new MultAll(ParType::MDOUBLE));

		std::vector<double> oneVec(sampleSize, 1.0);
		std::shared_ptr<AbsParameter> one(new MultiDouble("one",oneVec));
		newTree->createHead("Amplitude"+suffix, one);
		std::cout<<newTree->head()->to_str(10)<<std::endl;
		return newTree;
	}
};

#endif
