/*
 * Resonance.hpp
 *
 *  Created on: Mar 3, 2016
 *      Author: weidenka
 */

#ifndef CORE_RESONANCE_HPP_
#define CORE_RESONANCE_HPP_

#include <vector>
#include <memory>

#include <boost/iterator/filter_iterator.hpp>

#include "Core/AbsParameter.hpp"
#include "Core/Parameter.hpp"
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"

class Resonance
{
public:
	Resonance() {};
	virtual ~Resonance() {};

	//! Clone function
	virtual Resonance* Clone(std::string newName="") const = 0;

	//! Get resonance name
	virtual std::string GetName() const = 0;

	//! Set resonance name
	virtual void SetName(std::string name) = 0;

	//! Get phase of resonance name
	virtual bool GetEnable() const = 0;

	//! Implementation of interface for streaming info about the strategy
	virtual std::string to_str() const = 0;

	//!Get id of variable A
	virtual unsigned int GetVarIdA() const = 0;

	//!Get id of variable B
	virtual unsigned int GetVarIdB() const = 0;

	//!Set id of variable A
	virtual void SetVarIdA(unsigned int id) = 0;

	//!Set id of variable B
	virtual void SetVarIdB(unsigned int id) = 0;

	//! Get total integral for resonance
	virtual double GetIntegral() const = 0;

	//! Get total integral for resonance
	virtual double GetTotalIntegral() const = 0;

	//! Get coefficient
	virtual std::complex<double> GetCoefficient() const = 0;

	//! Get magnitude of resonance name
	virtual double GetMagnitude() const = 0;

	//! Get magnitude of resonance id
	virtual std::shared_ptr<DoubleParameter> GetMagnitudePar() = 0;

	//! Get phase of resonance name
	virtual double GetPhase() const = 0 ;

	//! Get phase of resonance id
	virtual std::shared_ptr<DoubleParameter> GetPhasePar() = 0;

	//! Get resonance mass
	virtual double GetMass() const = 0;

	//! Get resonance mass
	virtual std::shared_ptr<DoubleParameter> GetMassPar() = 0;

	//! Get resonance mass
	virtual double GetWidth() const = 0;

	//! Get resonance spin
	virtual double GetSpin() const = 0;

	//! value of resonance at \param point
	virtual std::complex<double> Evaluate(dataPoint& point) = 0;

	virtual std::shared_ptr<FunctionTree> SetupTree(
			ParameterList& sample, ParameterList& toySample,std::string suffix) = 0;
};

struct resIsEnabled {
	bool operator()(std::shared_ptr<Resonance> r) { return r->GetEnable(); }
};
typedef boost::filter_iterator<resIsEnabled,
		std::vector<std::shared_ptr<Resonance> >::iterator> resonanceItr;

#endif /* CORE_RESONANCE_HPP_ */
