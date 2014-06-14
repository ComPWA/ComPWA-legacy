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

#ifndef _MINUITRESULT_HPP_
#define _MINUITRESULT_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
//#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Core/ParameterList.hpp"
#include "Core/TableFormater.hpp"
#include "Core/PhysConst.hpp"
#include "Core/FitResult.hpp"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FunctionMinimum.h"

using namespace ROOT::Minuit2;

class MinuitResult : public FitResult
{
public:
	MinuitResult() {};
	MinuitResult(FunctionMinimum result) { init(result); }
	void setResult(FunctionMinimum result){ init(result); }
	void setInitialLH(double iniLH){ initialLH = iniLH; }
	void setAmplitude(std::shared_ptr<Amplitude> newAmp);
	operator double() const { return finalLH; };
	double getResult(){return finalLH;}
	void fractions(std::ostream& out);

private:
	bool isValid; //result valid
	bool covPosDef; //covariance matrix pos.-def.
	bool hasValidParameters; //valid parameters
	bool hasValidCov; //valid covariance
	bool hasAccCov; //accurate covariance
	bool hasReachedCallLimit; //call limit reached
	bool hesseFailed; //hesse failed

	double errorDef;
	unsigned int nFcn;
	double initialLH;
	double finalLH;
	double exitCode;
	double edm; //estimated distance to minimum
	boost::numeric::ublas::symmetric_matrix<double,boost::numeric::ublas::upper> cov;
	boost::numeric::ublas::symmetric_matrix<double,boost::numeric::ublas::upper> corr;
	boost::numeric::ublas::matrix<double> fracError;
	std::vector<double> variance;
	std::vector<double> globalCC;
	void genOutput(std::ostream& out,std::string opt="");
	void genSimpleOutput(std::ostream& out);
	void init(FunctionMinimum);
};

#endif
