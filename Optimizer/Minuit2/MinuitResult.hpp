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
#include <boost/numeric/ublas/io.hpp>

#include "Core/ParameterList.hpp"
#include "Core/TableFormater.hpp"
#include "Core/PhysConst.hpp"
#include "Core/FitResult.hpp"
#include "Core/Logging.hpp"
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
	//! Set amplitude
	void setAmplitude(std::shared_ptr<Amplitude> a) { _amp = a; }
	//! Enable correct error estimation for fit fractions. Very time consuming!
	void setUseCorrelatedErrors(bool s, int nSets=200);
	//! Use tree for calculation of fit fractions
	void setUseTree(bool s);
	//! Set calculation of interference terms
	void setCalcInterference(bool b) { calcInterference = b; }
	//! Write list of fit parameters and list of fitfractions to XML file @filename
	virtual void writeXML(std::string filename);
	//! Write fit parameters, fit fractions and cov matrix as TeX to file @filename
	virtual void writeTeX(std::string filename);
	//! Any errors during minimization?
	virtual bool hasFailed();
	//! Initialize result with Minuit2::FunctionMinimum
	void init(FunctionMinimum);

protected:
	//! Calculate interference terms
	bool calcInterference;
	//! Should we calcualte fit fraction errors accurately?
	bool useCorrelatedErrors;
	//! calculate fractions using tree (if available)
	bool useTree;
	//! number of resonances in amplitude
	unsigned int nRes;

	//! Number of floating parameters
	int nFreeParameter;
	//! Number of events
	int nEvents;
	//====== MINUIT FIT RESULT =======
	bool isValid; //result valid
	bool covPosDef; //covariance matrix pos.-def.
	bool hasValidParameters; //valid parameters
	bool hasValidCov; //valid covariance
	bool hasAccCov; //accurate covariance
	bool hasReachedCallLimit; //call limit reached
	bool edmAboveMax;
	bool hesseFailed; //hesse failed
	double errorDef;
	unsigned int nFcn;
	double initialLH;
	double finalLH;
	double penalty;
	double edm; //estimated distance to minimum
	//! Covariance matrix
	std::vector<std::vector<double> > cov;
	//! Correlation matrix
	std::vector<std::vector<double> > corr;
	//! Global correlation coefficients
	std::vector<double> globalCC;

	//====== OUTPUT =====
	//! Simplified fit result output
	void genSimpleOutput(std::ostream& out);
	//! Full fit result output
	void genOutput(std::ostream& out,std::string opt="");
	//! Table with correlation matrix
	void printCorrelationMatrix(TableFormater* fracTable);
	//! Table with covariance matrix
	void printCovarianceMatrix(TableFormater* fracTable);

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
	virtual void calcFractionError();
	//! Number of parameter sets that are used to propagate the cov matrix through the normalization
	unsigned int correlatedErrors_numberOfSets;
	//! Smear ParameterList with a multidimensional gaussian and the cov matrix from the fit
	void smearParameterList(ParameterList&);
	//! Calculate information criterion AIC
	double calcAIC();
	//! Calculate information criterion BIC
	double calcBIC();
private:
	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		using namespace boost::serialization;
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FitResult);
		ar & BOOST_SERIALIZATION_NVP(calcInterference);
		ar & BOOST_SERIALIZATION_NVP(useCorrelatedErrors);
		ar & BOOST_SERIALIZATION_NVP(useTree);
		ar & BOOST_SERIALIZATION_NVP(nRes);
		ar & BOOST_SERIALIZATION_NVP(isValid);
		ar & BOOST_SERIALIZATION_NVP(covPosDef);
		ar & BOOST_SERIALIZATION_NVP(hasValidParameters);
		ar & BOOST_SERIALIZATION_NVP(hasValidCov);
		ar & BOOST_SERIALIZATION_NVP(hasAccCov);
		ar & BOOST_SERIALIZATION_NVP(hasReachedCallLimit);
		ar & BOOST_SERIALIZATION_NVP(edmAboveMax);
		ar & BOOST_SERIALIZATION_NVP(hesseFailed);
		ar & BOOST_SERIALIZATION_NVP(errorDef);
		ar & BOOST_SERIALIZATION_NVP(nFcn);
		ar & BOOST_SERIALIZATION_NVP(initialLH);
		ar & BOOST_SERIALIZATION_NVP(finalLH);
		ar & BOOST_SERIALIZATION_NVP(penalty);
		ar & BOOST_SERIALIZATION_NVP(nEvents);
		ar & BOOST_SERIALIZATION_NVP(edm);
		ar & BOOST_SERIALIZATION_NVP(cov);
		ar & BOOST_SERIALIZATION_NVP(corr);
		ar & BOOST_SERIALIZATION_NVP(globalCC);
		ar & BOOST_SERIALIZATION_NVP(nFreeParameter);
	}
};


#endif
