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
#include "Estimator/Estimator.hpp"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FunctionMinimum.h"

#include <gsl/gsl_rng.h>

using namespace ROOT::Minuit2;

class MinuitResult : public FitResult
{
public:
	MinuitResult() {};
	MinuitResult(std::shared_ptr<ControlParameter> esti, FunctionMinimum result);
	void setResult(std::shared_ptr<ControlParameter> esti,FunctionMinimum result);
	void setInitialLH(double iniLH){ initialLH = iniLH; }
	//! Convert to double and return final LH values
	operator double() const { return finalLH; };
	//! Return final likelihood value
	double getResult(){return finalLH;}
	//!Generate output stream with fit result table
	void fractions(std::ostream& out);
	//! Enable correct error estimation for fit fractions. Very time consuming!
	void setUseCorrelatedErrors(bool s) { useCorrelatedErrors = s; }

private:
	bool isValid; //result valid
	bool covPosDef; //covariance matrix pos.-def.
	bool hasValidParameters; //valid parameters
	bool hasValidCov; //valid covariance
	bool hasAccCov; //accurate covariance
	bool hasReachedCallLimit; //call limit reached
	bool hesseFailed; //hesse failed
	//! number of resonances in amplitude
	unsigned int nRes;

	gsl_rng* r;//! GSL random generator, used for multivariate gauss
	bool useCorrelatedErrors;

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
	void genSimpleOutput(std::ostream& out);
	void smearParameterList(ParameterList&);
	void genOutput(std::ostream& out,std::string opt="");
	void init(FunctionMinimum);
	std::shared_ptr<Estimator> estimator;

	/** Calculate fit fractions.
	 * Fractions are calculated using the formular:
	 * \f[
	 *  f_i = \frac{|c_i|^2 \int A_i A_i^*}{\int \sum c_l c_m^* A_l A_m}
	 * \f]
	 * The \f$c_i\f$ complex coefficienct of the amplitude and the denominatior is the integral over
	 * the whole amplitude.
	 *
	 * @param frac result with fit fractions for the single resonances
	 */
	void calcFraction(std::vector<double>& frac);
	/** Calculate errors on fit result
	 * Set @param assumeUnCorrelatedErrors to assume that the error of the fit parameter only depends
	 * on the error of the magnitude. The error of normalization due the the fit error on magnitudes
	 * and phases is ignored.
	 * If we want to calculate the errors correctly we have to generate a set of fit parameters that
	 * are smeard by a multidimensional gaussian and the covariance matrix of the fit. For every set
	 * we calculate the fit frations and calculate its mean. The can be a very time consuming method,
	 * especially if the function tree is not used.
	 *
	 * @param fracError result with errors
	 */
	void calcFractionError(std::vector<double>& fracError);
};

#endif
