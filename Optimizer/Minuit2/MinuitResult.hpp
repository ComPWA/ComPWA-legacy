// Copyright (c) 2014, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// MinuitResult. Implementation of FitResult.
///

#ifndef _MINUITRESULT_HPP_
#define _MINUITRESULT_HPP_

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <boost/serialization/export.hpp>

#include <Minuit2/FunctionMinimum.h>

#include "Core/FitResult.hpp"
#include "Core/Logging.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Properties.hpp"
#include "Core/TableFormater.hpp"
#include "Estimator/Estimator.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

class MinuitResult : public FitResult {
public:
  /// Default constructor needed for boost::serialization.
  MinuitResult();

  /// Constructor with estimator and result
  MinuitResult(std::shared_ptr<ComPWA::Estimator::Estimator> esti,
               ROOT::Minuit2::FunctionMinimum result);

  /// Set Minuit2 function minimum
  void setResult(std::shared_ptr<ComPWA::Estimator::Estimator> esti,
                 ROOT::Minuit2::FunctionMinimum result);

  /// Return final likelihood value
  double result() { return FinalLH; }

  /// Set initial likelihood value
  virtual void setInitialLH(double iniLH) { InitialLH = iniLH; }

  /// Get initial likelihood value
  virtual double initialLH() { return InitialLH; }

  /// Set final likelihood value
  virtual void setFinalLH(double iniLH) { FinalLH = iniLH; }

  /// Get final likelihood value
  virtual double finalLH() { return FinalLH; }

  /// Set true likelihood value
  virtual void setTrueLH(double iniLH) { TrueLH = iniLH; }

  /// Get true likelihood value
  virtual double trueLH() { return TrueLH; }

  /// Convert to double and return final LH values
  operator double() const { return FinalLH; }

  /// Set calculation of interference terms
  void setCalcInterference(bool b) { CalcInterference = b; }

  /// Get calculation of interference terms
  bool calcInterference() { return CalcInterference; }

  /// Write list of fit parameters and list of fitfractions to XML file @param
  /// filename
  virtual void writeXML(std::string filename);

  /// Write fit parameters, fit fractions and cov matrix as TeX to file @param
  /// filename
  virtual void writeTeX(std::string filename);

  /// Is minimum valid?
  virtual bool isValid() { return IsValid; }

  /// Number of free parameters
  virtual int ndf() { return NumFreeParameter; }

  /// Get covariance matrix
  virtual std::vector<std::vector<double>> covarianceMatrix() { return Cov; }

  /// Get correlation matrix
  virtual std::vector<std::vector<double>> correlationMatrix() { return Corr; }

  /// Get global correlation coefficiencts
  virtual std::vector<double> gobalCC() { return GlobalCC; }

  /// Get estimated distrance to minimum
  virtual double edm() { return Edm; }

protected:
  /// Initialize result with Minuit2::FunctionMinimum
  void init(ROOT::Minuit2::FunctionMinimum);

  /// Pointer to estimator
  std::shared_ptr<ComPWA::Estimator::Estimator> est;
  
  /// Calculate interference terms
  bool CalcInterference;

  /// Number of floating parameters
  unsigned int NumFreeParameter;

  //====== MINUIT FIT RESULT =======
  double corr(unsigned int n, unsigned int t) {
    std::cout << "WARNING: not sure if row and column are choose correctly!"
              << std::endl;
    if (n < Corr.size() && t < Corr.at(1).size() && t >= n)
      return Corr.at(n).at(t);
    else
      return -9000;
  };

  bool IsValid;             // result valid
  bool CovPosDef;           // covariance matrix pos.-def.
  bool HasValidParameters;  // valid parameters
  bool HasValidCov;         // valid covariance
  bool HasAccCov;           // accurate covariance
  bool HasReachedCallLimit; // call limit reached
  bool EdmAboveMax;
  bool HesseFailed; // hesse failed
  double ErrorDef;
  unsigned int NFcn;
  double InitialLH;
  double FinalLH;
  double TrueLH;
  double Edm; // estimated distance to minimum

  /// Covariance matrix
  std::vector<std::vector<double>> Cov;

  /// Correlation matrix
  std::vector<std::vector<double>> Corr;

  /// Global correlation coefficients
  std::vector<double> GlobalCC;

  /// Full fit result output
  void genOutput(std::ostream &out, std::string opt = "");

  /// Create table with interference terms for each AmpIntensity
  void createInterferenceTable(std::ostream &out,
                               std::shared_ptr<Intensity> amp);

  /// Table with correlation matrix
  void printCorrelationMatrix(TableFormater *fracTable);

  /// Table with covariance matrix
  void printCovarianceMatrix(TableFormater *fracTable);

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    using namespace boost::serialization;
    ar &BOOST_SERIALIZATION_BASE_OBJECT_NVP(FitResult);
    ar &BOOST_SERIALIZATION_NVP(CalcInterference);
    ar &BOOST_SERIALIZATION_NVP(IsValid);
    ar &BOOST_SERIALIZATION_NVP(CovPosDef);
    ar &BOOST_SERIALIZATION_NVP(HasValidParameters);
    ar &BOOST_SERIALIZATION_NVP(HasValidCov);
    ar &BOOST_SERIALIZATION_NVP(HasAccCov);
    ar &BOOST_SERIALIZATION_NVP(HasReachedCallLimit);
    ar &BOOST_SERIALIZATION_NVP(EdmAboveMax);
    ar &BOOST_SERIALIZATION_NVP(HesseFailed);
    ar &BOOST_SERIALIZATION_NVP(ErrorDef);
    ar &BOOST_SERIALIZATION_NVP(NFcn);
    ar &BOOST_SERIALIZATION_NVP(InitialLH);
    ar &BOOST_SERIALIZATION_NVP(FinalLH);
    ar &BOOST_SERIALIZATION_NVP(TrueLH);
    ar &BOOST_SERIALIZATION_NVP(Edm);
    ar &BOOST_SERIALIZATION_NVP(Cov);
    ar &BOOST_SERIALIZATION_NVP(Corr);
    ar &BOOST_SERIALIZATION_NVP(GlobalCC);
    ar &BOOST_SERIALIZATION_NVP(NumFreeParameter);
  }
};

} // namespace Minuit2
} // namespace Optimizer
} // namespace ComPWA

#endif
