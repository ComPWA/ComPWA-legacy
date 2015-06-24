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
using namespace boost::log;

class FitResult
{
public:
	FitResult():time(0){};
	virtual ~FitResult() {};
	//! Set amplitude
	virtual void setAmplitude(std::shared_ptr<Amplitude> a) { _amp = a; }
	//! Set list of initial parameters
	virtual void setInitialParameters(ParameterList iniPars){ initialParameters=iniPars; }
	//! Set list of final fit parameters
	virtual void setFinalParameters(ParameterList finPars){ finalParameters=finPars; }
	//! Set list of true parameters
	virtual void setTrueParameters(ParameterList truePars){ trueParameters=truePars; }
	//! Set value of likelihood with initial parameter
	virtual void setInitialLH(double iniLH){ }
	//! Get list of initial parameters
	virtual ParameterList getInitialParameters(){ return initialParameters; }
	//! Get list of final fit parameters
	virtual ParameterList getFinalParameters(){ return finalParameters; }
	//! Get list of true parameters
	virtual ParameterList getTrueParameters(){ return trueParameters; }
	//! Set processing time for minimization
	virtual void setTime(double t){ time = t; }
	//! Get processing time for minimization
	virtual double getTime(){ return time; }
	//! Get fit result (e.g. likelihood or chi2)
	virtual double getResult() =0;
	//! Table with fit parameters
	virtual void printFitParameters(TableFormater* tableResult);
	//! Table with fit fractions
	virtual void printFitFractions(TableFormater* fracTable);
	//! Getter function for fractions list. Make sure that fractions are calculated beforehand.
	virtual ParameterList& GetFractions() {	return fractionList; }

	//output
	virtual void print(std::string opt=""){
		std::stringstream s;
		genOutput(s,opt);
		std::string str = s.str();
		BOOST_LOG_TRIVIAL(info) << str;
	};
	virtual void writeTeX(std::string filename) {};
	virtual void writeXML(std::string filename) {};
	virtual void writeText(std::string filename) ;
	virtual void writeSimpleText(std::string filename) ;
	virtual operator double() const =0;
	friend std::ostream& operator<< (std::ostream &out, FitResult &fitres){ out<<fitres.getResult(); return out;};
	//! Any errors during minimization?
	virtual bool hasFailed(){ return 0; };

protected:
	virtual double shiftAngle(double v);
	virtual void genOutput(std::ostream& out,std::string opt="") = 0;
	virtual void genSimpleOutput(std::ostream& out);

	double time;

	ParameterList initialParameters;
	ParameterList finalParameters;
	ParameterList trueParameters;
	std::shared_ptr<Amplitude> _amp;

	//! Calculate fit fractions and its errors.
	virtual void calcFraction();
	/** Calculate fit fractions.
	 * Fractions are calculated using the formular:
	 * \f[
	 *  f_i = \frac{|c_i|^2 \int A_i A_i^*}{\int \sum c_l c_m^* A_l A_m}
	 * \f]
	 * The \f$c_i\f$ complex coefficienct of the amplitude and the denominatior is the integral over
	 * the whole amplitude.
	 *
	 * @param parList result with fit fractions for the single resonances
	 */
	virtual void calcFraction(ParameterList& parList);
	//! Calculate errors on fit result
	virtual void calcFractionError() {};
	//! List with fit fractions and errors
	ParameterList fractionList;
};


#endif
