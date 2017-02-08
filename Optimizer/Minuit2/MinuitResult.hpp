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
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

namespace ComPWA {
namespace Optimizer {
namespace Minuit2 {

class MinuitResult : public FitResult {
public:
  //! Default constructor
  MinuitResult();

  //! Constructor with estimator and result
  MinuitResult(std::shared_ptr<ControlParameter> esti,
               ROOT::Minuit2::FunctionMinimum result);

  //! Set Minuit2 function minimum
  void setResult(std::shared_ptr<ControlParameter> esti,
                 ROOT::Minuit2::FunctionMinimum result);

  //! Return final likelihood value
  double getResult() { return finalLH; }

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

  //! Set list of true parameters
  virtual void setTrueParameters(ParameterList truePars);

  //! Set list of initial parameters
  virtual void setInitialParameters(ParameterList iniPars);

  //! Convert to double and return final LH values
  operator double() const { return finalLH; }

  //! Set calculation of interference terms
  void SetCalcInterference(bool b) { calcInterference = b; }

  //! Get calculation of interference terms
  bool GetCalcInterference() { return calcInterference; }

  //! Write list of fit parameters and list of fitfractions to XML file @param
  //! filename
  virtual void writeXML(std::string filename);

  //! Write fit parameters, fit fractions and cov matrix as TeX to file @param
  //! filename
  virtual void writeTeX(std::string filename);

  //! Any errors during minimization?
  virtual bool hasFailed();

  //! Is minimum valid?
  virtual bool MinimumIsValid() { return isValid; }

  //! Number of free parameters
  virtual int GetNDF() { return nFreeParameter; }

  //! Get number of events
  virtual int GetNEvents() { return nEvents; }

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

  //! Get AIC
  virtual double GetAIC() { return AIC; }

  //! Get BIC
  virtual double GetBIC() { return BIC; }

  //! Get penalty scale
  virtual double GetPenaltyScale() { return penaltyScale; }

  //! Get penalty term
  virtual double GetPenalty() { return penalty; }

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
  std::shared_ptr<ComPWA::Estimator::Estimator> est;

  //====== MINUIT FIT RESULT =======
  double getCorr(unsigned int n, unsigned int t) {
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
  double penalty;
  double penaltyScale;
  double edm; // estimated distance to minimum
  //! Covariance matrix
  std::vector<std::vector<double>> cov;
  //! Correlation matrix
  std::vector<std::vector<double>> corr;
  //! Global correlation coefficients
  std::vector<double> globalCC;

  double AIC;
  double BIC;

  //====== OUTPUT =====
  //! Simplified fit result output
  void genSimpleOutput(std::ostream &out);

  //! Full fit result output
  void genOutput(std::ostream &out, std::string opt = "");

  //! Create table with interference terms for each amplitude
  void createInterferenceTable(std::ostream &out,
                               std::shared_ptr<Amplitude> amp);

  //! Table with correlation matrix
  void printCorrelationMatrix(TableFormater *fracTable);

  //! Table with covariance matrix
  void printCovarianceMatrix(TableFormater *fracTable);

  /** Calculate errors on fit result
   * The error of normalization due the the fit error on magnitudes and phases
   * is ignored.
   * If we want to calculate the errors correctly we have to generate a set of
   * fit parameters that
   * are smeard by a multidimensional gaussian and the covariance matrix of the
   * fit. For every set
   * we calculate the fit frations and calculate its mean. The can be a very
   * time consuming method,
   * especially if the function tree is not used.
   *
   * @param parList fit parameter
   * @param amp amplitude that was fitted
   * @param nSets number of sets used
   */
  virtual void calcFractionError(ParameterList &parList,
                                 std::shared_ptr<Amplitude> amp, int nSets);

  //! Calculate information criterion AIC
  double calcAIC(ParameterList &frac);

  //! Calculate information criterion BIC
  double calcBIC(ParameterList &frac);

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
    ar &BOOST_SERIALIZATION_NVP(penalty);
    ar &BOOST_SERIALIZATION_NVP(penaltyScale);
    ar &BOOST_SERIALIZATION_NVP(AIC);
    ar &BOOST_SERIALIZATION_NVP(BIC);
    ar &BOOST_SERIALIZATION_NVP(nEvents);
    ar &BOOST_SERIALIZATION_NVP(edm);
    ar &BOOST_SERIALIZATION_NVP(cov);
    ar &BOOST_SERIALIZATION_NVP(corr);
    ar &BOOST_SERIALIZATION_NVP(globalCC);
    ar &BOOST_SERIALIZATION_NVP(nFreeParameter);
  }
};

/************** HELPER FUNCTION FOR GSL VECTOR AND MATRIX *******************/
/** Print gsl_matrix **/
inline void gsl_matrix_print(const gsl_matrix *m) {
  for (size_t i = 0; i < m->size1; i++) {
    for (size_t j = 0; j < m->size2; j++) {
      printf("%g ", gsl_matrix_get(m, i, j));
    }
    printf("\n");
  }
};

/** Print gsl_vector **/
inline void gsl_vector_print(const gsl_vector *m) {
  for (size_t i = 0; i < m->size; i++) {
    std::printf("%g ", gsl_vector_get(m, i));
  }
  std::printf("\n");
};

/** Convert ParameterList to gsl vector **/
inline gsl_vector *gsl_parameterList2Vec(const ParameterList &list) {
  gsl_vector *tmp = gsl_vector_alloc(list.GetNDouble());
  unsigned int t = 0;
  for (unsigned int o = 0; o < list.GetNDouble(); o++) {
    std::shared_ptr<DoubleParameter> outPar = list.GetDoubleParameter(o);
    if (outPar->IsFixed())
      continue;
    gsl_vector_set(tmp, t, outPar->GetValue());
    t++;
  }
  // resize vector
  gsl_vector *vec = gsl_vector_alloc(t);
  for (unsigned int i = 0; i < vec->size; i++)
    gsl_vector_set(vec, i, gsl_vector_get(tmp, i));
  return vec;
};

/** Convert std::vector matrix to gsl matrix **/
inline gsl_matrix *
gsl_vecVec2Matrix(const std::vector<std::vector<double>> &m) {
  gsl_matrix *tmp = gsl_matrix_alloc(m.size(), m.at(0).size());
  for (size_t i = 0; i < tmp->size1; i++) {
    for (unsigned int j = 0; j < tmp->size2; j++) {
      gsl_matrix_set(tmp, i, j, m.at(i).at(j));
    }
  }
  return tmp;
};

/** Multivariate Gaussian using cholesky decomposition
 * A test application can be found at test/MultiVariateGaussianTestApp.cpp
 *
 * @param rnd Random number generator
 * @param vecSize Size of data vector
 * @param in Mean value(s)
 * @param cov Covariance matrix
 * @param res Resulting Vector
 */
inline void multivariateGaussian(const gsl_rng *rnd, const int vecSize,
                                 const gsl_vector *in, const gsl_matrix *cov,
                                 gsl_vector *res) {
  // Generate and fill temporary covariance matrix
  gsl_matrix *tmpM = gsl_matrix_alloc(vecSize, vecSize);
  gsl_matrix_memcpy(tmpM, cov);

  // Cholesky decomposition
  int status = gsl_linalg_cholesky_decomp(tmpM);
  if (status == GSL_EDOM)
    LOG(error) << "Decomposition has failed!";

  // Compute vector of random gaussian variables
  for (unsigned int i = 0; i < vecSize; i++)
    gsl_vector_set(res, i, gsl_ran_ugaussian(rnd));

  // Debug
  //	gsl_matrix_print(cov);
  //	gsl_vector_print(in);
  //	gsl_vector_print(res);
  //	gsl_matrix_print(tmpM);

  /*From the GNU GSL Documentation:
   * The function dtrmv compute the matrix-vector product x = op(A) x for the
   * triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans,
   * CblasTrans, CblasConjTrans. When Uplo is CblasUpper then the upper
   * triangle of A is used, and when Uplo is CblasLower then the lower
   * triangle of A is used. If Diag is CblasNonUnit then the diagonal of
   * the matrix is used, but if Diag is CblasUnit then the diagonal elements
   * of the matrix A are taken as unity and are not referenced.
   */
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, tmpM, res);

  gsl_vector_add(res, in);
  // free temporary object
  gsl_matrix_free(tmpM);
};
/************** HELPER FUNCTION FOR GSL VECTOR AND MATRIX *******************/

} /* namespace Minuit2 */
} /* namespace Optimizer */
} /* namespace ComPWA */

#endif
