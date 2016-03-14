//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//			Peter Weidenkaff
//-------------------------------------------------------------------------------
//! Optimizer Interface Base-Class.
/*! \class Optimizer
 * @file Optimizer.hpp
 * This class provides the interface to (external) optimization libraries or
 * routines. As it is pure virtual, one needs at least one implementation to
 * provide an optimizer for the analysis which varies free model-parameters. If
 * a new optimizer is derived from and fulfills this base-class, no change in
 * other modules are necessary to work with the new optimizer library or routine.
 */

#ifndef _FITRESULT_HPP_
#define _FITRESULT_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/log/trivial.hpp>

#include "Core/Amplitude.hpp"
#include "Core/ParameterList.hpp"
#include "Core/TableFormater.hpp"
#include "Core/PhysConst.hpp"

namespace ComPWA {

class FitResult
{
public:
	FitResult():time(0){};
	virtual ~FitResult() {};
	void setInitialParameters(ParameterList iniPars){ initialParameters=iniPars; }
	void setFinalParameters(ParameterList finPars){ finalParameters=finPars; }
	void setTrueParameters(ParameterList truePars){ trueParameters=truePars; }
	void setInitialLH(double iniLH){ }
	ParameterList getInitialParameters(){ return initialParameters; }
	ParameterList getFinalParameters(){ return finalParameters; }
	ParameterList getTrueParameters(){ return trueParameters; }
	void setTime(double t){ time = t; }
	double getTime(){ return time; }
	virtual double getResult() =0;

	//output
	virtual void print(std::string opt=""){
		std::stringstream s;
		genOutput(s,opt);
		std::string str = s.str();
		BOOST_LOG_TRIVIAL(info) << str;
	};
	//! Table with fit parameters
	virtual void printFitParameters(TableFormater* tableResult);
	virtual void writeTeX(std::string filename) {};
	virtual void writeXML(std::string filename) {};
	virtual void writeText(std::string filename) ;
	virtual void writeSimpleText(std::string filename) ;
	virtual operator double() const =0;
	friend std::ostream& operator<< (std::ostream &out, FitResult &fitres){ out<<fitres.getResult(); return out;};
	//! Any errors during minimization?
	virtual bool hasFailed(){};

protected:
	virtual double shiftAngle(double v);
	virtual void genOutput(std::ostream& out,std::string opt="") = 0;
	virtual void genSimpleOutput(std::ostream& out);

	double time;

	ParameterList initialParameters;
	ParameterList finalParameters;
	ParameterList trueParameters;
	std::shared_ptr<Amplitude> _amp;
};

} /* namespace ComPWA */

#endif
