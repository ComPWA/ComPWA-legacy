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
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/export.hpp>

#include "Core/Amplitude.hpp"
#include "Core/ParameterList.hpp"
#include "Core/TableFormater.hpp"
#include "Core/PhysConst.hpp"

namespace ComPWA {


class FitResult
{
public:
	FitResult() : time(0), nSetsFractionError(0) { };
	virtual ~FitResult() {};
	//! Set single amplitude. Assume that only one amplitude is used!
	virtual void SetAmplitude(std::shared_ptr<Amplitude> a) {
		_ampVec.clear();
		AddAmplitude(a);
	}
	//! Add amplitude
	virtual void AddAmplitude(std::shared_ptr<Amplitude> a) {
		_ampVec.push_back(a);
	}
	//! Set amplitude vector
	virtual void SetAmplitude( std::vector<std::shared_ptr<Amplitude> > vec) {
		_ampVec = vec;
	}
	//! Set fraction list
	virtual void setFractions(ParameterList ini){ fractionList=ini; }
	//! Set list of initial parameters
	virtual void setInitialParameters(ParameterList iniPars){ initialParameters.DeepCopy(iniPars); }
	//! Set list of final fit parameters
	virtual void setFinalParameters(ParameterList finPars);
	//! Set list of true parameters
	virtual void setTrueParameters(ParameterList truePars){ trueParameters.DeepCopy(truePars); }
	//! Set value of likelihood with initial parameter
	virtual void SetInitialLH(double iniLH){ }
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
	virtual double getCorr(unsigned int n, unsigned int t) {return -9000;};

	//! Table with fit parameters
	virtual void printFitParameters(TableFormater* tableResult);
	//! Table with fit fractions
	virtual void printFitFractions(TableFormater* tab);
	//! Table with fit fractions
	virtual void printFitFractions(TableFormater* tab,
			std::shared_ptr<Amplitude> amp, int nErrorSets=0);
	//! Getter function for fractions list. Make sure that fractions are calculated beforehand.
	virtual ParameterList& getFractions() {	return fractionList; }

	//! Enable correct error estimation for fit fractions. Very time consuming!
	void setUseCorrelatedErrors(int nSets=200);

	//! Print fit result
	virtual void print(std::string opt="");

	virtual void writeTeX(std::string filename) {};
	virtual void writeXML(std::string filename) {};
	virtual void writeText(std::string filename) ;
	virtual void writeSimpleText(std::string filename) ;
	virtual operator double() const =0;
	friend std::ostream& operator<< (std::ostream &out, FitResult &fitres){
		out<<fitres.getResult();
		return out;
	};
	//! Any errors during minimization?
	virtual bool hasFailed(){ return 0; };

protected:
	virtual double shiftAngle(double v);
	virtual void genOutput(std::ostream& out,std::string opt="") = 0;
	virtual void genSimpleOutput(std::ostream& out);

	//! Time for minimization
	double time;
	//! Initial list of parameters
	ParameterList initialParameters;
	//! Final list of parameters
	ParameterList finalParameters;
	//! True list of parameters
	ParameterList trueParameters;
	//! Fit amplitude (can't be serialized)
	std::vector<std::shared_ptr<Amplitude> > _ampVec;

	//! Number of parameter sets that are used to propagate the cov matrix through the normalization
	int nSetsFractionError;

	/** Calculate fit fractions.
	 * Fractions are calculated using the formular:
	 * \f[
	 *  f_i = \frac{|c_i|^2 \int A_i A_i^*}{\int \sum c_l c_m^* A_l A_m}
	 * \f]
	 * The \f$c_i\f$ complex coefficienct of the amplitude and the denominatior is the integral over
	 * the whole amplitude.
	 *
	 * @param parList result with fit fractions for the single resonances
	 * @param nSets Precise error calcucation using @nSets Monte-Carlo events
	 */
//	virtual void calcFraction(ParameterList& parList, int nSets=0);

	//! Calculate fit fractions and its errors.
//	static void calcFraction(
//			ParameterList& parList, std::shared_ptr<Amplitude> amp);

	//! Calculate errors on fit result
	virtual void calcFractionError(ParameterList& parList,
			std::shared_ptr<Amplitude> amp, int nSets=200) = 0;

	//! List with fit fractions and errors
	ParameterList fractionList;
	double sumFractions;
	double sumFractionsError;


private:
	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(time);
		ar & BOOST_SERIALIZATION_NVP(initialParameters);
		ar & BOOST_SERIALIZATION_NVP(finalParameters);
		ar & BOOST_SERIALIZATION_NVP(trueParameters);
		ar & BOOST_SERIALIZATION_NVP(fractionList);
		ar & BOOST_SERIALIZATION_NVP(sumFractions);
		ar & BOOST_SERIALIZATION_NVP(sumFractionsError);
		ar & BOOST_SERIALIZATION_NVP(nSetsFractionError);
	}
};
BOOST_SERIALIZATION_ASSUME_ABSTRACT( FitResult );
} /* namespace ComPWA */

#endif
