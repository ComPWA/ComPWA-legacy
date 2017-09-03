// Copyright (c) 2014, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Optimizer Interface Base-Class.
/*! \class Optimizer
 * @file Optimizer.hpp
 * This class provides the interface to (external) optimization libraries or
 * routines. As it is pure virtual, one needs at least one implementation to
 * provide an optimizer for the analysis which varies free model-parameters. If
 * a new optimizer is derived from and fulfills this base-class, no change in
 * other modules are necessary to work with the new optimizer library or
 * routine.
 */

#ifndef _MINUITRESULT_HPP_
#define _MINUITRESULT_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include <Minuit2/FunctionMinimum.h>

#include "Core/ParameterList.hpp"
#include "Core/TableFormater.hpp"
#include "Core/Properties.hpp"
#include "Core/FitResult.hpp"
#include "Core/Logging.hpp"
#include "Core/Estimator.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

class MinuitResult : public FitResult {
public:
  //! Default constructor
  MinuitResult();

  //! Constructor with estimator and result
  MinuitResult(std::shared_ptr<ComPWA::IEstimator> esti,
               ROOT::Minuit2::FunctionMinimum result);

  //! Set Minuit2 function minimum
  void setResult(std::shared_ptr<ComPWA::IEstimator> esti,
                 ROOT::Minuit2::FunctionMinimum result);

  //! Return final likelihood value
  double GetResult() { return finalLH; }

  //! Set initial likelihood value
  virtual void SetInitialLH(double iniLH) { initialLH = iniLH; }
  //! Get initial likelihood value
  virtual double GetInitialLH() { return initialLH; }
  //! Set final likelihood value
  virtual void SetFinalLH(double iniLH) { finalLH = iniLH; }
  //! Get final likelihood value
  virtual double GetFinalLH() { return finalLH; }
  //! Set true likelihood value
  virtual void SetTrueLH(double iniLH) { trueLH = iniLH; }
  //! Get true likelihood value
  virtual double GetTrueLH() { return trueLH; }

  //! Convert to double and return final LH values
  operator double() const { return finalLH; }

  //! Set calculation of interference terms
  void SetCalcInterference(bool b) { calcInterference = b; }

  //! Get calculation of interference terms
  bool GetCalcInterference() { return calcInterference; }

  //! Write list of fit parameters and list of fitfractions to XML file @param
  //! filename
  virtual void WriteXML(std::string filename);

  //! Write fit parameters, fit fractions and cov matrix as TeX to file @param
  //! filename
  virtual void WriteTeX(std::string filename);

  //! Any errors during minimization?
  virtual bool HasFailed();

  //! Is minimum valid?
  virtual bool MinimumIsValid() { return isValid; }

  //! Number of free parameters
  virtual int GetNDF() { return nFreeParameter; }

  //! Get covariance matrix
  virtual std::vector<std::vector<double>> GetCovarianceMatrix() { return cov; }
  
  //! Get correlation matrix
  virtual std::vector<std::vector<double>> GetCorrelationMatrix() {
    return corr;
  }
  
  //! Get global correlation coefficiencts
  virtual std::vector<double> GetGlobalCC() { return globalCC; }
  
  //! Get estimated distrance to minimum
  virtual double GetEDM() { return edm; }

protected:
  //! Initialize result with Minuit2::FunctionMinimum
  void init(ROOT::Minuit2::FunctionMinimum);

  //! Calculate interference terms
  bool calcInterference;

  //! Number of floating parameters
  int nFreeParameter;

  //! Number of events
  int nEvents;

  //! Pointer to estimator
  std::shared_ptr<ComPWA::IEstimator> est;

  //====== MINUIT FIT RESULT =======
  double GetCorr(unsigned int n, unsigned int t) {
    std::cout << "WARNING: not sure if row and column are choose correctly!"
              << std::endl;
    if (n < corr.size() && t < corr.at(1).size() && t >= n)
      return corr.at(n).at(t);
    else
      return -9000;
  };

private:
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

  //====== OUTPUT =====
  //! Simplified fit result output
  void genSimpleOutput(std::ostream &out);

  //! Full fit result output
  void genOutput(std::ostream &out, std::string opt = "");

  //! Create table with interference terms for each AmpIntensity
  void createInterferenceTable(std::ostream &out,
                               std::shared_ptr<AmpIntensity> amp);

  //! Table with correlation matrix
  void PrintCorrelationMatrix(TableFormater *fracTable);

  //! Table with covariance matrix
  void PrintCovarianceMatrix(TableFormater *fracTable);

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    using namespace boost::serialization;
    ar &BOOST_SERIALIZATION_BASE_OBJECT_NVP(FitResult);
    ar &BOOST_SERIALIZATION_NVP(calcInterference);
    ar &BOOST_SERIALIZATION_NVP(isValid);
    ar &BOOST_SERIALIZATION_NVP(covPosDef);
    ar &BOOST_SERIALIZATION_NVP(hasValidParameters);
    ar &BOOST_SERIALIZATION_NVP(hasValidCov);
    ar &BOOST_SERIALIZATION_NVP(hasAccCov);
    ar &BOOST_SERIALIZATION_NVP(hasReachedCallLimit);
    ar &BOOST_SERIALIZATION_NVP(edmAboveMax);
    ar &BOOST_SERIALIZATION_NVP(hesseFailed);
    ar &BOOST_SERIALIZATION_NVP(errorDef);
    ar &BOOST_SERIALIZATION_NVP(nFcn);
    ar &BOOST_SERIALIZATION_NVP(initialLH);
    ar &BOOST_SERIALIZATION_NVP(finalLH);
    ar &BOOST_SERIALIZATION_NVP(trueLH);
//    ar &BOOST_SERIALIZATION_NVP(penalty);
//    ar &BOOST_SERIALIZATION_NVP(penaltyScale);
//    ar &BOOST_SERIALIZATION_NVP(AIC);
//    ar &BOOST_SERIALIZATION_NVP(BIC);
    ar &BOOST_SERIALIZATION_NVP(nEvents);
    ar &BOOST_SERIALIZATION_NVP(edm);
    ar &BOOST_SERIALIZATION_NVP(cov);
    ar &BOOST_SERIALIZATION_NVP(corr);
    ar &BOOST_SERIALIZATION_NVP(globalCC);
    ar &BOOST_SERIALIZATION_NVP(nFreeParameter);
  }
};

} /* namespace Minuit2 */
} /* namespace Optimizer */
} /* namespace ComPWA */

#endif
