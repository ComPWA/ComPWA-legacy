// Some utility functions -*- C++ -*-
/* Copyright 2008 Mike Williams (mwill@jlab.org)
 *
 * This file is part of qft++.
 *
 * qft++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * qft++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with qft++.  If not, see <http://www.gnu.org/licenses/>.
 */
// Author: Mike Williams
#ifndef _WIGNER_D_HPP
#define _WIGNER_D_HPP
//_____________________________________________________________________________
/** @file
 *  @brief Some utility functions.
 */
//_____________________________________________________________________________
#include <cctype>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Core/Spin.hpp"

namespace ComPWA {
namespace QFT {

inline int factorial(int __i) {
  int f = 1;
  if ((__i == 0) || (__i == 1))
    f = 1;
  else {
    while (__i > 0) {
      f = f * __i;
      __i--;
    }
  }
  return f;
}

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)

/// Utility function used by Clebsch
inline double dfact(double __x) {
  if ((__x < 0.00001) && (__x >= 0.0))
    return 1.;
  if (__x < 0)
    return 0.;
  return __x * dfact(__x - 1.);
}

inline double Clebsch(const Spin &_j1, const Spin &_j2, const Spin &_J) {
  // convert to pure integers (each 2*spin)
  auto mag = _j1.getMagnitude();
  int j1 = 2 * mag.getNumerator() / mag.getDenominator();
  auto zproj = _j1.getProjection();
  int m1 = 2 * zproj.getNumerator() / zproj.getDenominator();
  mag = _j2.getMagnitude();
  int j2 = 2 * mag.getNumerator() / mag.getDenominator();
  zproj = _j2.getProjection();
  int m2 = 2 * zproj.getNumerator() / zproj.getDenominator();
  mag = _J.getMagnitude();
  int J = 2 * mag.getNumerator() / mag.getDenominator();
  zproj = _J.getProjection();
  int M = 2 * zproj.getNumerator() / zproj.getDenominator();

  if ((m1 + m2) != M)
    return 0.;
  if (m1 > j1)
    return 0;
  if (m2 > j2)
    return 0;

  double n0, n1, n2, n3, n4, n5, d0, d1, d2, d3, d4, A, exp;
  int nu = 0;

  double sum = 0;
  while (((d3 = (j1 - j2 - M) / 2 + nu) < 0) ||
         ((n2 = (j1 - m1) / 2 + nu) < 0)) {
    nu++;
  }
  while (((d1 = (J - j1 + j2) / 2 - nu) >= 0) &&
         ((d2 = (J + M) / 2 - nu) >= 0) &&
         ((n1 = (j2 + J + m1) / 2 - nu) >= 0)) {
    d3 = ((j1 - j2 - M) / 2 + nu);
    n2 = ((j1 - m1) / 2 + nu);
    d0 = dfact((double)nu);
    exp = nu + (j2 + m2) / 2;
    n0 = (double)pow(-1., exp);
    sum += ((n0 * dfact(n1) * dfact(n2)) /
            (d0 * dfact(d1) * dfact(d2) * dfact(d3)));
    nu++;
  }

  if (sum == 0)
    return 0;

  n0 = J + 1;
  n1 = dfact((double)(J + j1 - j2) / 2);
  n2 = dfact((double)(J - j1 + j2) / 2);
  n3 = dfact((double)(j1 + j2 - J) / 2);
  n4 = dfact((double)(J + M) / 2);
  n5 = dfact((J - M) / 2);

  d0 = dfact((double)(j1 + j2 + J) / 2 + 1);
  d1 = dfact((double)(j1 - m1) / 2);
  d2 = dfact((double)(j1 + m1) / 2);
  d3 = dfact((double)(j2 - m2) / 2);
  d4 = dfact((double)(j2 + m2) / 2);

  A = ((double)(n0 * n1 * n2 * n3 * n4 * n5)) /
      ((double)(d0 * d1 * d2 * d3 * d4));

  return sqrt(A) * sum;
}

//_____________________________________________________________________________
/** Returns \f$d^{j}_{m,n}(\beta)\f$.
 * Calculates the Wigner d-functions.
 *
 *  This function was copied from one written by D.P. Weygand. It shouldn't be
 *  used for spins higher than about 11/2 or 8.
 */
inline double Wigner_d(const Fraction &_J, const Fraction &_M,
                       const Fraction &_N, double _beta) {

  int J = 2 * _J.getNumerator() / _J.getDenominator();
  int M = 2 * _M.getNumerator() / _M.getDenominator();
  int N = 2 * _N.getNumerator() / _N.getDenominator();
  double beta = _beta;

  int temp_M, k, k_low, k_hi;
  double const_term = 0.0, sum_term = 0.0, d = 1.0;
  int m_p_n, j_p_m, j_p_n, j_m_m, j_m_n;
  int kmn1, kmn2, jmnk, jmk, jnk;
  double kk;

  if (J < 0 || abs(M) > J || abs(N) > J) {
    std::cerr << std::endl;
    std::cerr << "d: you have entered an illegal number for J, M, N."
              << std::endl;
    std::cerr
        << "Must follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J."
        << std::endl;
    std::cerr << "J = " << J << " M = " << M << " N = " << N << std::endl;
    return 0.;
  }

  if (beta < 0) {
    beta = fabs(beta);
    temp_M = M;
    M = N;
    N = temp_M;
  }

  m_p_n = (M + N) / 2;
  j_p_m = (J + M) / 2;
  j_m_m = (J - M) / 2;
  j_p_n = (J + N) / 2;
  j_m_n = (J - N) / 2;

  kk = (double)factorial(j_p_m) * (double)factorial(j_m_m) *
       (double)factorial(j_p_n) * (double)factorial(j_m_n);
  const_term = pow((-1.0), (j_p_m)) * sqrt(kk);

  k_low = MAX(0, m_p_n);
  k_hi = MIN(j_p_m, j_p_n);

  for (k = k_low; k <= k_hi; k++) {

    kmn1 = 2 * k - (M + N) / 2;
    jmnk = J + (M + N) / 2 - 2 * k;
    jmk = (J + M) / 2 - k;
    jnk = (J + N) / 2 - k;
    kmn2 = k - (M + N) / 2;

    sum_term +=
        pow((-1.0), (k)) *
        ((pow(cos(beta / 2.0), kmn1)) * (pow(sin(beta / 2.0), jmnk))) /
        (factorial(k) * factorial(jmk) * factorial(jnk) * factorial(kmn2));
  }

  d = const_term * sum_term;
  return d;
}

} // namespace QFT
} // namespace ComPWA
#endif
