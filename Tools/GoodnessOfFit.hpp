#ifndef COMPWA_TOOLS_GOF_HPP_
#define COMPWA_TOOLS_GOF_HPP_

#include <cmath>

/// Calculate Akaike Information Criterium AICc
/// https://en.wikipedia.org/wiki/Akaike_information_criterion
/// @param LH         assuming -log(L) is passed
/// @param sampleSize number of points in the data sample
/// @param nFreePar   number of free parameters
inline double calculateAIC(double LH, int sampleSize, int nFreePar) {
  // Classic AIC
  double aic = 2.0 * LH + 2.0 * nFreePar;
  // AICc is a modified version which includes the sample size
  aic += 2.0 * (nFreePar * nFreePar + nFreePar) / (sampleSize - nFreePar - 1);
  return aic;
}

/// Calculate Beyesian Information Criterium BIC
/// https://en.wikipedia.org/wiki/Bayesian_information_criterion
/// @param LH    assuming -log(L) is passed
/// @param sampleSize number of points in the data sample
/// @param nFreePar   number of free parameters
inline double calculateBIC(double LH, int sampleSize, int nFreePar) {
  return std::log(sampleSize) * nFreePar + 2.0 * LH;
}
#endif
