
//
//  FitFractions.hpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 06.05.17.
//
//
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

inline ComPWA::DoubleParameter
CalculateFitFraction(std::shared_ptr<ComPWA::AmpIntensity> intens,
                     std::shared_ptr<std::vector<dataPoint>> sample,
                     const std::pair<std::string, std::string> def) {

  std::shared_ptr<ComPWA::AmpIntensity> numer, denom;
  if (def.first == intens->Name())
    numer = intens;
  else
    numer = intens->GetComponent(def.first);

  if (def.second == intens->Name())
    denom = intens;
  else
    denom = intens->GetComponent(def.second);

  double integral_numerator = ComPWA::Tools::Integral(numer, sample);
  double integral_denominator = Tools::Integral(denom, sample);
  double ffVal = integral_numerator / integral_denominator;
  LOG(trace) << "CalculateFitFraction() | Result for (" << def.first << "/"
             << def.second << ") is " << integral_numerator << "/"
             << integral_denominator << "=" << ffVal;
  return DoubleParameter(def.first, ffVal, 0.0);
}

inline ComPWA::ParameterList
CalculateFitFractions(std::shared_ptr<ComPWA::AmpIntensity> intens,
                      std::shared_ptr<std::vector<dataPoint>> sample,
                      std::vector<std::pair<std::string, std::string>> defs) {
  ComPWA::ParameterList ffList;
  for (auto i : defs) {
    auto par = CalculateFitFraction(intens, sample, i);
    ffList.AddParameter(std::make_shared<ComPWA::DoubleParameter>(par));
  }
  return ffList;
}

/** Calculate fit fractions.
 * Fractions are calculated using the formular:
 * \f[
 *  f_i = \frac{|c_i|^2 \int A_i A_i^*}{\int \sum c_l c_m^* A_l A_m}
 * \f]
 * The \f$c_i\f$ complex coefficienct of the amplitude and the denominatior is
 * the integral over
 * the whole amplitude.
 *
 * @param parList result with fit fractions for the single resonances
 * @param amp amplitude with the resonances which share the fit fractions
 * @param nSets Precise error calcucation using @param nSets Monte-Carlo
 * events
 */
/** Calculate errors on fit fractions
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
 * @param amp AmpIntensity that was fitted
 * @param nSets number of sets used
 */
inline void
CalcFractionError(ParameterList &parameters,
                  std::vector<std::vector<double>> covariance,
                  ParameterList &ffList, std::shared_ptr<AmpIntensity> intens,
                  std::shared_ptr<std::vector<dataPoint>> sample, int nSets,
                  std::vector<std::pair<std::string, std::string>> defs) {
  if (nSets <= 0)
    return;
  if (!parameters.GetNDouble())
    return;
  LOG(info) << "Calculating errors of fit fractions using " << nSets
            << " sets of parameters...";

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

  std::vector<ParameterList> fracVect;
  progressBar bar(nSets);
  std::stringstream outFraction;
  //	for(unsigned int i=0; i<nSets; i++){
  int i = 0;
  while (i < nSets) {
    bool error = 0;
    bar.nextEvent();
    gsl_vector *gslNewPar = gsl_vector_alloc(nFreeParameter);
    // generate set of smeared parameters
    multivariateGaussian(rnd, nFreeParameter, gslFinalPar, gslCov, gslNewPar);
    gsl_vector_print(gslNewPar);

    // deep copy of finalParameters
    ParameterList newPar;
    newPar.DeepCopy(parameters);

    std::size_t t = 0;
    for (std::size_t o = 0; o < newPar.GetNDouble(); o++) {
      std::shared_ptr<DoubleParameter> outPar = newPar.GetDoubleParameter(o);
      if (outPar->IsFixed())
        continue;
      // set floating values to smeared values
      try { // catch out-of-bound
        outPar->SetValue(gslNewPar->data[t]);
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
      //      Amplitude::UpdateAmpParameterList(_ampVec, newPar);
    } catch (ParameterOutOfBound &ex) {
      continue;
    }
    fracVect.push_back(CalculateFitFractions(intens, sample, defs));
    i++;

    /******* DEBUGGING *******/
    //			if(i==0){
    //				for(int t=0; t<newPar.GetNDouble();
    // t++){
    //					if(
    // newPar.GetDoubleParameter(t)->IsFixed())
    // continue;
    //					outFraction <<
    // newPar.GetDoubleParameter(t)->GetName()<<":";
    //				}
    //				for(int t=0; t<tmp.GetNDouble(); t++)
    //					outFraction <<
    // tmp.GetDoubleParameter(t)->GetName()<<":";
    //				outFraction << "norm" << std::endl;
    //			}
    //			for(int t=0; t<newPar.GetNDouble(); t++){
    //				if(
    // newPar.GetDoubleParameter(t)->IsFixed())
    // continue;
    //				outFraction <<
    // newPar.GetDoubleParameter(t)->GetValue()<<"
    //";
    //			}
    //			for(int t=0; t<tmp.GetNDouble(); t++)
    //				outFraction <<
    // tmp.GetDoubleParameter(t)->GetValue()<<"
    //";
    //			double norm = _amp->GetIntegral();
    //			outFraction << norm;
    //			outFraction << std::endl;
    /******* DEBUGGING *******/
  }
  LOG(info) << " ------- " << outFraction.str();

  // free objects
  gsl_vector_free(gslFinalPar);
  gsl_matrix_free(gslCov);
  gsl_rng_free(rnd);

  int nRes = parameters.GetNDouble();
  // Calculate standard deviation
  for (unsigned int o = 0; o < nRes; o++) {
    double mean = 0, sqSum = 0., stdev = 0;
    for (unsigned int i = 0; i < fracVect.size(); i++) {
      double tmp = fracVect.at(i).GetDoubleParameter(o)->GetValue();
      mean += tmp;
      sqSum += tmp * tmp;
    }
    unsigned int s = fracVect.size();
    sqSum /= s;
    mean /= s;
    // this is cross-checked with the RMS of the distribution
    stdev = std::sqrt(sqSum - mean * mean);
    parameters.GetDoubleParameter(o)->SetError(stdev);
  }

  // Set correct fit result
  //  intens->UpdateParameters(parameters);
  return;
}
} /* namespace Tools */
} /* namespace ComPWA */
#endif /* FitFractions_h */
