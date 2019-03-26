// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Optimizer Interface Base-Class.
/*! \class GenevaResult
 * @file GenevaResult.hpp
 * This class contains FitResults from Geneva Optimizations
 */

#ifndef _GENEVARESULT_HPP_
#define _GENEVARESULT_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Core/ParameterList.hpp"
#include "Core/TableFormater.hpp"
#include "Core/FitResult.hpp"

#include "Optimizer/Geneva/GStartIndividual.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Geneva {

class GenevaResult : public FitResult
{
public:
	GenevaResult() : FitResult() {};
	void setResult(std::shared_ptr<Gem::Geneva::GStartIndividual> result){ init(result); }
	operator double() const { return finalLH; };
	double result(){ return finalLH; }

protected:
	  bool isValid;             // result valid
	  bool covPosDef;           // covariance matrix pos.-def.
	  bool hasValidParameters;  // valid parameters
	  bool hasValidCov;         // valid covariance
	  bool hasAccCov;           // accurate covariance
	  bool hasReachedCallLimit; // call limit reached
	  bool edmAboveMax;
	  bool hesseFailed; // hesse failed
	  double errorDef;
	  unsigned int nFcn;
	  double initialLH;
	  double finalLH;
	  double trueLH;
	//  double penalty;
	//  double penaltyScale;
	  double edm; // estimated distance to minimum
	  //! Covariance matrix
	  std::vector<std::vector<double>> cov;
	  //! Correlation matrix
	  std::vector<std::vector<double>> corr;
	  //! Global correlation coefficients
	  std::vector<double> globalCC;

	  //! Table with correlation matrix
	  void printCorrelationMatrix(TableFormater *fracTable);

	  //! Table with covariance matrix
	  void printCovarianceMatrix(TableFormater *fracTable);

	  void genOutput(std::ostream& out, std::string opt="");

	  void init(std::shared_ptr<Gem::Geneva::GStartIndividual> min);

	  //virtual void calcFractionError(ParameterList& parList,
	  //		std::shared_ptr<Amplitude> amp, int nSets=200) { };
};

} /* namespace Geneva */
} /* namespace Optimizer */
} /* namespace ComPWA */

#endif
