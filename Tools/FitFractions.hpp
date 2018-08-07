// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <memory>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "Core/ProgressBar.hpp"
#include "Tools/Integration.hpp"

#ifndef FitFractions_h
#define FitFractions_h

namespace ComPWA {
namespace Tools {

/// Print gsl_matrix
inline void gsl_matrix_print(const gsl_matrix *m) {
  for (size_t i = 0; i < m->size1; i++) {
    for (size_t j = 0; j < m->size2; j++) {
      printf("%g ", gsl_matrix_get(m, i, j));
    }
    printf("\n");
  }
};

/// Print gsl_vector
inline void gsl_vector_print(const gsl_vector *m) {
  for (size_t i = 0; i < m->size; i++) {
    std::printf("%g ", gsl_vector_get(m, i));
  }
  std::printf("\n");
};

/// Convert ParameterList to gsl vector
inline gsl_vector *gsl_parameterList2Vec(const ParameterList &list) {
  size_t n = list.doubleParameters().size();
  gsl_vector *tmp = gsl_vector_alloc(n);
  size_t t = 0;
  for (size_t o = 0; o < n; o++) {
    auto outPar = list.doubleParameter(o);
    if (outPar->isFixed())
      continue;
    gsl_vector_set(tmp, t, outPar->value());
    t++;
  }
  // resize vector
  gsl_vector *vec = gsl_vector_alloc(t);
  for (unsigned int i = 0; i < vec->size; i++)
    gsl_vector_set(vec, i, gsl_vector_get(tmp, i));
  return vec;
};

/// Convert std::vector matrix to gsl matrix
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

/// Multivariate Gaussian using cholesky decomposition
/// A test application can be found at test/MultiVariateGaussianTestApp.cpp
///
/// \param rnd Random number generator
/// \param vecSize Size of data vector
/// \param in Mean value(s)
/// \param cov Covariance matrix
/// \param res Resulting Vector
inline void multivariateGaussian(const gsl_rng *rnd, const unsigned int vecSize,
                                 const gsl_vector *in, const gsl_matrix *cov,
                                 gsl_vector *res) {
  // Generate and fill temporary covariance matrix
  gsl_matrix *tmpM = gsl_matrix_alloc(vecSize, vecSize);
  gsl_matrix_memcpy(tmpM, cov);

  // Cholesky decomposition
  int status = gsl_linalg_cholesky_decomp(tmpM);
  if (status == GSL_EDOM)
    LOG(ERROR) << "Decomposition has failed!";

  // Compute vector of random gaussian variables
  for (unsigned int i = 0; i < vecSize; i++)
    gsl_vector_set(res, i, gsl_ran_ugaussian(rnd));

  // Debug
  //  gsl_matrix_print(cov);
  //  gsl_vector_print(in);
  //  gsl_vector_print(res);
  //  gsl_matrix_print(tmpM);

  // From the GNU GSL Documentation:
  // The function dtrmv compute the matrix-vector product x = op(A) x for the
  // triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans,
  // CblasTrans, CblasConjTrans. When Uplo is CblasUpper then the upper
  // triangle of A is used, and when Uplo is CblasLower then the lower
  // triangle of A is used. If Diag is CblasNonUnit then the diagonal of
  // the matrix is used, but if Diag is CblasUnit then the diagonal elements
  // of the matrix A are taken as unity and are not referenced.
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, tmpM, res);

  gsl_vector_add(res, in);
  // free temporary object
  gsl_matrix_free(tmpM);
};

inline ComPWA::FitParameter
CalculateFitFraction(std::shared_ptr<ComPWA::Kinematics> kin,
                     std::shared_ptr<ComPWA::AmpIntensity> intens,
                     std::shared_ptr<std::vector<DataPoint>> sample,
                     const std::pair<std::string, std::string> def) {

  double phspVolume = kin->phspVolume();
  
  std::shared_ptr<ComPWA::AmpIntensity> denom;
  if (def.second == intens->name())
    denom = intens;
  else
    denom = intens->component(def.second);

  double integral_denominator =
      ComPWA::Tools::Integral(denom, sample, phspVolume);
  
  std::shared_ptr<ComPWA::AmpIntensity> numer;
  if (def.first == denom->name())
    numer = denom;
  else
    numer = denom->component(def.first);
  double integral_numerator =
      ComPWA::Tools::Integral(numer, sample, phspVolume);

  double ffVal = integral_numerator / integral_denominator;
  LOG(TRACE) << "CalculateFitFraction() | Result for (" << def.first << "/"
             << def.second << ") is " << integral_numerator << "/"
             << integral_denominator << "=" << ffVal;
  return FitParameter(def.first, ffVal, 0.0);
}

/// Calculate fit fractions.
/// Fractions are calculated using the formular:
/// \f[
///  f_i = \frac{|c_i|^2 \int A_i A_i^*}{\int \sum c_l c_m^* A_l A_m}
/// \f]
/// The \f$c_i\f$ complex coefficienct of the amplitude and the denominatior is
/// the integral over
/// the whole amplitude.
inline ComPWA::ParameterList
CalculateFitFractions(std::shared_ptr<ComPWA::Kinematics> kin,
                      std::shared_ptr<ComPWA::AmpIntensity> intens,
                      std::shared_ptr<std::vector<DataPoint>> sample,
                      std::vector<std::pair<std::string, std::string>> defs) {
  ComPWA::ParameterList ffList;
  for (auto i : defs) {
    auto par = CalculateFitFraction(kin, intens, sample, i);
    ffList.addParameter(std::make_shared<ComPWA::FitParameter>(par));
  }
  return ffList;
}

/// Calculate errors on fit fractions.
/// The error of normalization due the the fit error on magnitudes and phases
/// is ignored. If we want to calculate the errors correctly we have to
/// generate a set of fit parameters that are smeard by a multidimensional
/// gaussian and the covariance matrix of the fit. For every set we calculate
/// the fit frations and calculate its mean. The can be a very time consuming
/// method, especially if the FunctionTree is not used.
inline void CalcFractionError(
    ParameterList &parameters, std::vector<std::vector<double>> covariance,
    ParameterList &ffList, std::shared_ptr<ComPWA::Kinematics> kin,
    std::shared_ptr<AmpIntensity> intens,
    std::shared_ptr<std::vector<DataPoint>> sample, int nSets,
    std::vector<std::pair<std::string, std::string>> defs) {
  LOG(INFO)
      << "CalcFractionError() | Calculating errors of fit fractions using "
      << nSets << " sets of parameters...";

  if (nSets <= 0)
    return;
  if (!parameters.doubleParameters().size())
    return;
  if (parameters.doubleParameters().size() != covariance.size())
      return;
  
  // Check wheather covariance matrix is set to zero.
  bool leave = true;
  for(unsigned int i=0; i<covariance.size(); ++i)
    for (unsigned int j=0; j < covariance.size(); ++j)
      if (covariance.at(i).at(j) != 0.)
        leave = false;
  if (leave)
    LOG(ERROR) << "CalcFractionError() | Covariance matrix is zero "
                  "(everywhere)! We skip further "
                  "calculation of fit fraction errors.";

  // Setting up random number generator
  const gsl_rng_type *T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  gsl_rng *rnd = gsl_rng_alloc(T);

  // convert to GSL objects
  gsl_vector *gslFinalPar = gsl_parameterList2Vec(parameters);
  gsl_matrix *gslCov = gsl_vecVec2Matrix(covariance);
  gsl_matrix_print(gslCov); // DEBUG

  int nFreeParameter = covariance.size();

  ParameterList originalPar;
  originalPar.DeepCopy(parameters);
  
  std::vector<ParameterList> fracVect;
  ProgressBar bar(nSets);
  int i = 0;
  while (i < nSets) {
    bool error = 0;
    gsl_vector *gslNewPar = gsl_vector_alloc(nFreeParameter);
    // generate set of smeared parameters
    multivariateGaussian(rnd, nFreeParameter, gslFinalPar, gslCov, gslNewPar);
    //gsl_vector_print(gslNewPar); // Debug

    // deep copy of finalParameters
    ParameterList newPar;
    newPar.DeepCopy(parameters);

    std::size_t t = 0;
    for (std::size_t o = 0; o < newPar.doubleParameters().size(); o++) {
      std::shared_ptr<FitParameter> outPar = newPar.doubleParameter(o);
      if (outPar->isFixed())
        continue;
      // set floating values to smeared values
      try { // catch out-of-bound
        outPar->setValue(gslNewPar->data[t]);
      } catch (ParameterOutOfBound &ex) {
        error = 1;
      }
      t++;
    }
    if (error)
      continue; // skip this set if one parameter is out of bound

    // free vector
    gsl_vector_free(gslNewPar);
    // update amplitude with smeared parameters
    try {
      intens->updateParameters(newPar);
    } catch (ParameterOutOfBound &ex) {
      continue;
    }
    fracVect.push_back(CalculateFitFractions(kin, intens, sample, defs));
    for( auto i : defs)
      std::cout << i.first << " " << i.second << std::endl;
    bar.next();
    i++;

    /******* DEBUGGING *******/
    //      if(i==0){
    //        for(int t=0; t<newPar.GetNDouble();
    // t++){
    //          if(
    // newPar.GetFitParameter(t)->IsFixed())
    // continue;
    //          outFraction <<
    // newPar.GetFitParameter(t)->name()<<":";
    //        }
    //        for(int t=0; t<tmp.GetNDouble(); t++)
    //          outFraction <<
    // tmp.GetFitParameter(t)->name()<<":";
    //        outFraction << "norm" << std::endl;
    //      }
    //      for(int t=0; t<newPar.GetNDouble(); t++){
    //        if(
    // newPar.GetFitParameter(t)->IsFixed())
    // continue;
    //        outFraction <<
    // newPar.GetFitParameter(t)->value()<<"
    //";
    //      }
    //      for(int t=0; t<tmp.GetNDouble(); t++)
    //        outFraction <<
    // tmp.GetFitParameter(t)->value()<<"
    //";
    //      double norm = _amp->GetIntegral();
    //      outFraction << norm;
    //      outFraction << std::endl;
    /******* DEBUGGING *******/
  }

  // free objects
  gsl_vector_free(gslFinalPar);
  gsl_matrix_free(gslCov);
  gsl_rng_free(rnd);

  // Calculate standard deviation
  for (unsigned int o = 0; o < ffList.doubleParameters().size(); ++o) {
    double mean = 0, sqSum = 0., stdev = 0;
    for (unsigned int i = 0; i < fracVect.size(); ++i) {
      double tmp = fracVect.at(i).doubleParameter(o)->value();
      mean += tmp;
      sqSum += tmp * tmp;
    }
    unsigned int s = fracVect.size();
    sqSum /= s;
    mean /= s;
    // Equivalent to RMS of the distribution
    stdev = std::sqrt(sqSum - mean * mean);
    ffList.doubleParameter(o)->setError(stdev);
  }

  // Set correct fit result
  intens->updateParameters(originalPar);
  return;
}

} // ns::Tools
} // ns::ComPWA

#endif
