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

#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "DataReader/Data.hpp"

#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Generator.hpp"
class Amplitude
{

public:

	Amplitude(){}

	virtual ~Amplitude()
	{ /* nothing */	}

	virtual Amplitude* Clone() = 0;

	//! set efficiency
	virtual void setEfficiency(std::shared_ptr<Efficiency> eff) {};
	virtual const double integral() =0;
	virtual const double integral(ParameterList& par) =0;
	virtual const double normalization() =0;
	virtual const double normalization(ParameterList& par) =0;
	virtual double getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen) = 0;
	virtual double getMaxVal(std::shared_ptr<Generator> gen) = 0;

	virtual const ParameterList& intensity(dataPoint& point, ParameterList& par) =0;
	virtual const ParameterList& intensity(dataPoint& point) =0;
	virtual const ParameterList& intensityNoEff(dataPoint& point) =0;
	virtual const ParameterList& intensity(std::vector<double> point, ParameterList& par) =0;

	virtual void setParameterList(ParameterList& par) =0;
	virtual bool copyParameterList(ParameterList& par) =0;

	virtual void printAmps() = 0;
	virtual void printFractions() = 0;

	/** Integral value of amplitude in certain boundary
	 * Used for plotting a projection of a function in \p var1 in bin [\p min1, \p min2]. In this
	 * case we have to integrate over an other variable \p var2
	 * @param var1 first variable
	 * @param min1 minimal value of first variable
	 * @param max1 maximal value of first variable
	 * @param var2 second variable
	 * @param min2 minimal value of second variable
	 * @param max2 maximal value of second variable
	 * @return
	 */
	virtual double getIntValue(std::string var1, double min1, double max1, std::string var2, double min2, double max2) = 0;

	//---------- get resonance parameters -------------
	//! Number of resonances
	virtual unsigned int getNumberOfResonances() { return 0; }
	//! convert resonance \param name to id
	virtual int getIdOfResonance(std::string name){ return 0;}
	//! convert resonance \param id to name
	virtual std::string getNameOfResonance(unsigned int id){ return std::string("muh");}
	//! get total integral for resonance \param id
	virtual double getAmpIntegral(unsigned int id) { return -999; };
	//! get total integral for resonance \param name
	virtual double getAmpIntegral(std::string name) { return -999; };
	//! Get magnitude of resonance name
	virtual double getAmpMagnitude(std::string name) {return -999;};
	//! Get magnitude of resonance id
	virtual double getAmpMagnitude(unsigned int id) {return -999;};
	//! Get phase of resonance name
	virtual double getAmpPhase(std::string name) {return -999;};
	//! Get phase of resonance id
	virtual double getAmpPhase(unsigned int id) {return -999;};

	//---------- related to FunctionTree -------------
	//! Check of tree is available
	virtual bool hasTree(){ return 0; }
	//! Getter function for basic amp tree
	virtual std::shared_ptr<FunctionTree> getAmpTree(allMasses&,allMasses&, std::string){
		return std::shared_ptr<FunctionTree>();
	}

protected:
	ParameterList result;


};

class GaussAmp : public Amplitude
{
public:
	GaussAmp(const char *name, DoubleParameter _resMass, DoubleParameter _resWidth){
		params.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter(_resMass)));
		params.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter(_resWidth)));
		initialise();
	}

	GaussAmp(const char *name, double _resMass, double _resWidth){
		params.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("mass",_resMass)));
		params.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("width",_resWidth)));
		initialise();
	}
	virtual void initialise() {
		result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("GaussAmpResult")));
		if(Kinematics::instance()->getVarNames().size()!=1)
			throw std::runtime_error("GaussAmp::initialize() | this amplitude is for two body decays only!");
	};
	//! Clone function
	virtual GaussAmp* Clone(){
		return (new GaussAmp(*this));
	}
	virtual bool copyParameterList(ParameterList& outPar){
		outPar = ParameterList(params);
		return true;
	}

	virtual double getIntValue(std::string var1, double min1, double max1,
			std::string var2, double min2, double max2) { return 0; }
	/**! Integral from -inf to inf
	 * @return
	 */
	virtual const double integral(){
		return (params.GetDoubleParameter(1)->GetValue() * std::sqrt(2*M_PI));
	}
	virtual const double integral(ParameterList& par){
		setParameterList(par);
		return integral();
	}
	virtual const double normalization() { return integral(); }
	virtual const double normalization(ParameterList& par){
		setParameterList(par);
		return normalization();
	}
	virtual double getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen){
		setParameterList(par);
		return getMaxVal(gen);
	}
	virtual double getMaxVal(std::shared_ptr<Generator> gen){
		double mass = params.GetDoubleParameter(0)->GetValue();
		std::vector<double> m; m.push_back(mass*mass);
		dataPoint p(m);
		intensity(p);
		return (result.GetParameterValue(0));
	}
	virtual const ParameterList& intensity(dataPoint& point, ParameterList& par){
		setParameterList(par);
		return intensity(point);
	}
	virtual const ParameterList& intensity(dataPoint& point) {
		double mass = params.GetDoubleParameter(0)->GetValue();
		double width = params.GetDoubleParameter(1)->GetValue();
		double sqrtS = std::sqrt(point.getVal(0));

		std::complex<double> gaus(std::exp(-1*(sqrtS-mass)*(sqrtS-mass)/width/width/2.),0);
		if(gaus.real() != gaus.real())
			BOOST_LOG_TRIVIAL(error)<<"GaussAmp::intensity() | result NaN!";
		result.SetParameterValue(0,std::norm(gaus));
		return result;
	}
	virtual const ParameterList& intensityNoEff(dataPoint& point){ return intensity(point); }
	virtual const ParameterList& intensity(std::vector<double> point, ParameterList& par){
		setParameterList(par);
		dataPoint dataP(point);
		return intensity(dataP);
	}
	void setParameterList(ParameterList& par){
		//parameters varied by Minimization algorithm
		if(par.GetNDouble()!=params.GetNDouble())
			throw std::runtime_error("setParameterList(): size of parameter lists don't match");
		//Should we compared the parameter names? String comparison is slow
		for(unsigned int i=0; i<params.GetNDouble(); i++)
			params.GetDoubleParameter(i)->UpdateParameter(par.GetDoubleParameter(i));
		return;
	}

	virtual void printAmps() { };
	virtual void printFractions() { };

protected:
	std::string _name;
	ParameterList params;
	//	std::shared_ptr<DoubleParameter> _mR;
	//	std::shared_ptr<DoubleParameter> _resWidth;


};

/**! UnitAmp
 *
 * Example implementation of Amplitude with the function value 1.0 at all points in PHSP. It is used
 * to test the likelihood normalization.
 */
class UnitAmp : public Amplitude
{
public:
	UnitAmp() : eff_(std::shared_ptr<Efficiency>(new UnitEfficiency())){
		result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("AmpSumResult")));
	}

	virtual ~UnitAmp()	{ /* nothing */	}

	virtual Amplitude* Clone() {;
	return (new UnitAmp(*this));
	}

	//! set efficiency
	virtual void setEfficiency(std::shared_ptr<Efficiency> eff) { eff_ = eff; }
	virtual const double integral() {
		return Kinematics::instance()->getPhspVolume();
	}
	virtual const double integral(ParameterList& par) { return integral(); }
	virtual const double normalization() {
		BOOST_LOG_TRIVIAL(info) << "UnitAmp::normalization() | normalization not implemented!";
		return 1;
	}
	virtual const double normalization(ParameterList& par) { return normalization(); }
	virtual double getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen) { return 1; }
	virtual double getMaxVal(std::shared_ptr<Generator> gen) { return 1; }

	virtual const ParameterList& intensity(dataPoint& point, ParameterList& par) {
		return intensity(point);
	}
	virtual const ParameterList& intensity(dataPoint& point) {
		result.SetParameterValue(0,eff_->evaluate(point));
		return result;
	}
	virtual const ParameterList& intensityNoEff(dataPoint& point) {
		result.SetParameterValue(0,1.0);
		return result;
	}
	virtual const ParameterList& intensity(std::vector<double> point, ParameterList& par){
		dataPoint dataP(point);
		return intensity(dataP);
	}

	virtual void setParameterList(ParameterList& par) { }
	virtual bool copyParameterList(ParameterList& par) { return 1; }

	virtual void printAmps() { }
	virtual void printFractions() { }

	virtual double getIntValue(std::string var1, double min1, double max1,
			std::string var2, double min2, double max2) { return 0; }

	//---------- related to FunctionTree -------------
	//! Check of tree is available
	virtual bool hasTree(){ return 1; }
	//! Getter function for basic amp tree
	virtual std::shared_ptr<FunctionTree> getAmpTree(allMasses& theMasses,
			allMasses& toyPhspSample, std::string suffix){
		return setupBasicTree(theMasses,toyPhspSample, suffix);
	}

protected:
	std::shared_ptr<Efficiency> eff_;

	/**Setup Basic Tree
	 *
	 * @param theMasses data sample
	 * @param toyPhspSample sample of flat toy MC events for normalization of the resonances
	 * @param opt Which tree should be created? "data" data Tree, "norm" normalization tree
	 * with efficiency corrected toy phsp sample or "normAcc" normalization tree with sample
	 * of accepted flat phsp events
	 */
	std::shared_ptr<FunctionTree> setupBasicTree(allMasses& theMasses,
			allMasses& toyPhspSample, std::string suffix="") {
		BOOST_LOG_TRIVIAL(debug) << "UnitAmp::setupBasicTree() generating new tree!";
		if(theMasses.nEvents==0){
			BOOST_LOG_TRIVIAL(error) << "UnitAmp::setupBasicTree() data sample empty!";
			return std::shared_ptr<FunctionTree>();
		}
		std::shared_ptr<FunctionTree> newTree(new FunctionTree());
		//std::shared_ptr<MultAll> mmultDStrat(new MultAll(ParType::MDOUBLE));

		std::vector<double> oneVec(theMasses.nEvents, 1.0);
		std::shared_ptr<MultiDouble> one(new MultiDouble("one",oneVec));
		newTree->createHead("Amplitude"+suffix, one);
//		newTree->createHead("Amplitude"+suffix, 3.0);
//		if(!newTree->head()){
//			std::cout<<"asfdasfasdfas"<<std::endl;
//			exit(1);
//		}
		std::cout<<newTree->head()->to_str(10)<<std::endl;
		return newTree;
	}
};

#endif
