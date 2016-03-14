//-------------------------------------------------------------------------------
// Copyright (c) 2013 michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Data Interface Base-Class.
/*! \class Data
 * @file Data.hpp
 * This class provides the interface to experimental data. As it is pure virtual,
 * one needs at least one implementation to provide data for the other modules. If
 * a new reader is derived from and fulfills this base-class, no change in other
 * modules are necessary to work with the new dataset.
 */

#ifndef DATA_HPP_
#define DATA_HPP_

#include <vector>
#include <string>
#include <memory>

#include "Core/Event.hpp"
#include "Core/Generator.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {
namespace DataReader {

class Data
{

public:

	Data(){ }
	virtual ~Data()	{ /* nothing */	}

	virtual void Add(Data& otherSample) = 0;
	virtual void pushEvent(const Event&) =0;

	virtual const unsigned int getNEvents() const =0;
	virtual Event& getEvent(const int) =0;
	virtual std::vector<Event> getEvents() = 0;
	virtual allMasses getMasses(const unsigned int startEvent=0, unsigned int nEvents=0) = 0;

	virtual const unsigned int getNBins() const =0;
	virtual const int getBin(const int, double&, double&) =0; //TODO: BinDataTyp, dynamic dimension

	//! Remove all events outside PHSP
	virtual void reduceToPhsp() {};
	//! Select only first @param newSize events from full sample
	virtual void reduce(unsigned int newSize) = 0;
	//! Select random subset of events
	virtual std::shared_ptr<Data> rndSubSet(unsigned int size, std::shared_ptr<Generator> gen) = 0;

	//! Set efficiency value for all stored events. Efficiency is taken from Efficiency object.
	virtual void setEfficiency(std::shared_ptr<Efficiency> eff) { };
	//! Reset effciencies of all events
	virtual void resetEfficiency(double e=1.) { };

	virtual void resetWeights(double w=1.) = 0;
	virtual bool hasWeights() = 0;
	virtual void Clear() = 0;

	virtual void writeData(std::string file="", std::string trName="") =0;
};

} /* namespace DataReader */
} /* namespace ComPWA */

#endif
