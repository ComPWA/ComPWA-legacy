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

#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "DataReader/Data.hpp"

#include "Core/DataPoint.hpp"
#include "Core/Generator.hpp"
class Amplitude
{

public:

	Amplitude(){}

	virtual ~Amplitude()
	{ /* nothing */	}

	virtual const double integral() =0;
	virtual const double integral(ParameterList& par) =0;
	virtual double getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen) = 0;
	virtual double getMaxVal(std::shared_ptr<Generator> gen) = 0;
	//virtual const double volume() =0;

	virtual const ParameterList& intensity(dataPoint& point, ParameterList& par) =0;
	virtual const ParameterList& intensity(dataPoint& point) =0;
	virtual const ParameterList& intensityNoEff(dataPoint& point) =0;
	virtual const ParameterList& intensity(std::vector<double> point, ParameterList& par) =0;

	virtual const bool fillStartParVec(ParameterList& outPar) =0;
	virtual void setParameterList(ParameterList& par) =0;
	virtual void copyParameterList(ParameterList& par) {}

	virtual void printAmps() = 0;
	virtual void printFractions() = 0;
	virtual unsigned int getNumberOfResonances() { return 0; }

	//! get total integral for resonance \param id
	virtual double getTotalIntegral(unsigned int id) { return -999; };
	//! get total integral for resonance \param name
	virtual double getTotalIntegral(std::string name) { return -999; };
	//! convert resonance \param name to id
	virtual int getIdOfResonance(std::string name){ return 0;}
	//! convert resonance \param id to name
	virtual std::string getNameOfResonance(unsigned int id){ return std::string("muh");}
	virtual double getMagnitude(std::string name) {return -999;};
	virtual double getMagnitude(unsigned int id) {return -999;};
	virtual double getPhase(std::string name) {return -999;};
	virtual double getPhase(unsigned int id) {return -999;};
	virtual double getSpin(std::string name) {return -999;};
	virtual double getSpin(unsigned int id) {return -999;};
	virtual double getFraction(std::string name) = 0;
	virtual double getFraction(unsigned int id) = 0;
	virtual double getIntValue(std::string var1, double min1, double max1, std::string var2, double min2, double max2) = 0;
	virtual Amplitude* Clone() = 0;

	//! Check of tree is available
	virtual bool hasTree(){ return 0; }
	//! Getter function for basic amp tree
	virtual std::shared_ptr<FunctionTree> getAmpTree(allMasses&,allMasses&, std::string){
		return std::shared_ptr<FunctionTree>();
	}


	/* OBSOLETE SECTION ONLY FOR TESTING */
	virtual std::shared_ptr<FunctionTree> functionTree(allMasses& theMasses, allMasses& toyPhspSample) {
		//if not implemented, return NULL-pointer
		return std::shared_ptr<FunctionTree>();
	}
	virtual void resetTree() {
		//if not implemented, return NULL-pointer
		return;
	}
	virtual std::shared_ptr<FunctionTree> phspTree(allMasses& accPhspSample, allMasses& toyPhspSample) {
		//if not implemented, return NULL-pointer
		return std::shared_ptr<FunctionTree>();
	}
	virtual std::shared_ptr<FunctionTree> phspTree(allMasses& toyPhspSample) {
		//if not implemented, return NULL-pointer
		return std::shared_ptr<FunctionTree>();
	}
	//! Getter function for function tree
	virtual std::shared_ptr<FunctionTree> getTree(){ return std::shared_ptr<FunctionTree>(); }
	//! Getter function for phsp tree
	virtual std::shared_ptr<FunctionTree> getPhspTree(){ return std::shared_ptr<FunctionTree>(); }
	/* OBSOLETE SECTION ONLY FOR TESTING */

protected:
	ParameterList result;


};

#endif
