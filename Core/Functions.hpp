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
//! Functions to be used in FuntionTree.
/*! \class Strategy
 * \class AddAll
 * \class MultAll
 * \class PowerTwo
 * @file Functions.hpp
 * This file contains Functions implementing the Strategy interface so they
 * can be used inside a node of the FuntionTree to calculate the node-value.
 * In addition to the simple functions provided here, the interface can also
 * be used at other places to provide functions for the FunctionTree.
 */

#ifndef _FUNCTIONS_HPP_
#define _FUNCTIONS_HPP_

#include <vector>
#include <complex>
#include <math.h>

#include "Core/Exceptions.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"
#include "Core/AbsParameter.hpp"
#include "Core/DataPoint.hpp"

class Strategy
{
public:
	//! Constructor
	Strategy(ParType in):checkType(in){
	};
	virtual ~Strategy(){}

	//! friend function to stream parameter information to output
	/*!
	 * Declaring the stream-operator << as friend allows to stream parameter
	 * information to the output as easily as a generic type.
	 * \sa make_str(), to_str()
	 */
	friend std::ostream& operator<<( std::ostream& out,
			std::shared_ptr<Strategy> b ){
		return out << b->to_str();
	}

	//! friend function to stream parameter information to output
	/*!
	 * Declaring the stream-operator << as friend allows to stream parameter
	 * information to the output as easily as a generic type.
	 * \sa make_str(), to_str()
	 */
	friend std::ostream& operator<<( std::ostream& out, const Strategy& b ){
		return out << b.to_str();
	}

	//! Get ParType
	virtual const ParType OutType() const {
		return checkType;
	}

	//! Pure Virtual interface for streaming info about the strategy
	virtual const std::string to_str() const =0;

	//! Pure Virtual interface for executing a strategy
	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out) = 0;

protected:
	ParType checkType;
};

class Inverse : public Strategy
{
public:
	Inverse(ParType in):Strategy(in){
	};
	virtual ~Inverse() {};

	virtual const std::string to_str() const{
		return "+";
	}

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out);
};

class SquareRoot : public Strategy
{
public:
	SquareRoot(ParType in):Strategy(in){
	};
	virtual ~SquareRoot(){}

	virtual const std::string to_str() const{
		return "+";
	}

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out);
};

class AddAll : public Strategy
{
public:
	AddAll(ParType in):Strategy(in){
	};
	virtual ~AddAll(){}

	virtual const std::string to_str() const{
		return "+";
	}

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out);
};

class MultAll : public Strategy
{
public:
	MultAll(ParType in):Strategy(in){
	};
	virtual ~MultAll(){}

	virtual const std::string to_str() const{
		return "*";
	};

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out);
};

class LogOf : public Strategy
{
public:
	LogOf(ParType in):Strategy(in){
	};
	virtual ~LogOf(){};

	virtual const std::string to_str() const{
		return "Log";
	};

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out);
};

class Complexify : public Strategy
{
public:
	Complexify(ParType in):Strategy(in){
	};
	virtual ~Complexify(){}

	virtual const std::string to_str() const{
		return "MakeComplex";
	};

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out);
};

class AbsSquare : public Strategy
{
public:
	AbsSquare(ParType in):Strategy(in){
	};
	virtual ~AbsSquare(){}

	virtual const std::string to_str() const{
		return "||^2";
	};

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out);
};

class Power : public Strategy
{
public:
	Power(ParType in):Strategy(in){
	};
	virtual ~Power(){}

	virtual const std::string to_str() const{
		return "^";
	}

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out);
};

#endif
