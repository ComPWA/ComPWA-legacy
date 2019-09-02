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
#include <cmath>
#include <iostream>
#include <complex>
#include <string>
#include <cctype>
#include <vector>
#include <map>

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

inline double Clebsch(const double _j1, const double _m1, const double _j2,
               const double _m2, const double _J, const double _M) {
  // convert to pure integers (each 2*spin)
  int j1 = (int)(2. * _j1);
  int m1 = (int)(2. * _m1);
  int j2 = (int)(2. * _j2);
  int m2 = (int)(2. * _m2);
  int J = (int)(2. * _J);
  int M = (int)(2. * _M);

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
inline double Wigner_d(const double _J, const double _M, const double _N,
                double _beta) {

  int J = 2 * _J;
  int M = 2 * _M;
  int N = 2 * _N;
  double beta = _beta;

  int temp_M, k, k_low, k_hi;
  double const_term = 0.0, sum_term = 0.0, d = 1.0;
  int m_p_n, j_p_m, j_p_n, j_m_m, j_m_n;
  int kmn1, kmn2, jmnk, jmk, jnk;
  double kk;

  if (J < 0 || abs(M) > J || abs(N) > J) {
    std::cerr << std::endl;
    std::cerr << "d: you have entered an illegal number for J, M, N." << std::endl;
    std::cerr << "Must follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J."
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

//_____________________________________________________________________________
/** Returns the Wigner D-function \f$D^{j}_{m,n}(\alpha,\beta,\gamma)\f$.
 *  This function uses the single spin Wigner_d function, thus it shares its
 *  limitations regarding higher spins.
 */
inline std::complex<double> Wigner_D(double __alpha, double __beta, double __gamma,
                                const double __j, const double __m,
                                const double __n) {
  std::complex<double> i(0., 1.);
  return exp(-i * __m * __alpha + __n * __gamma) *
         Wigner_d(__j, __m, __n, __beta);
}

//_____________________________________________________________________________
/** Fills @a d w/ \f$d^{j}_{m,n}(\beta)\f$ for all spins from 0 to @a jmax
 *  (integer spin) or 1/2 to @a jmax (half-integer spin). Calculations are
 *  performed recursively, thus these are valid to much higher spin than the
 *  single-spin version.
 */
/*inline void Wigner_d(const Spin &__jmax, double __beta,
              std::map<Spin, std::map<Spin, std::map<Spin, double>>> &__d) {
  __d.clear();
  Spin jmin, one_half = 1 / 2.;
  if ((double)__jmax == 0.) {
    __d[0][0][0] = 1.;
    return;
  }
  double cb = cos(__beta);
  double sb = sin(__beta);
  // j=1 d's
  std::map<Spin, std::map<Spin, double>> d1;
  d1[1][1] = (1 + cb) / 2.;
  d1[1][0] = -sb / sqrt(2.);
  d1[1][-1] = (1 - cb) / 2.;
  d1[0][1] = sb / sqrt(2.);
  d1[0][0] = cb;
  d1[0][-1] = -sb / sqrt(2.);
  d1[-1][1] = (1 - cb) / 2.;
  d1[-1][0] = sb / sqrt(2.);
  d1[-1][-1] = (1 + cb) / 2.;

  if (__jmax.GetDenominator() == 1) { // integral spins
    __d[0][0][0] = 1.0;
    if ((double)__jmax == 0.)
      return;
    __d[1][1][1] = d1[1][1];
    __d[1][1][0] = d1[1][0];
    __d[1][1][-1] = d1[1][-1];
    __d[1][0][1] = d1[0][1];
    __d[1][0][0] = d1[0][0];
    __d[1][0][-1] = d1[0][-1];
    __d[1][-1][1] = d1[-1][1];
    __d[1][-1][0] = d1[-1][0];
    __d[1][-1][-1] = d1[-1][-1];
    if ((double)__jmax == 1.)
      return;
    jmin = 2.;
  } else { // half-integral spins
    __d[one_half][one_half][one_half] = cos(__beta / 2.);
    __d[one_half][one_half][(-1) * (int)one_half] = -sin(__beta / 2.);
    __d[one_half][(-1) * (int)one_half][one_half] = sin(__beta / 2.);
    __d[one_half][(-1) * (int)one_half][(-1) * (int)one_half] =
        cos(__beta / 2.);
    if (__jmax == one_half)
      return;
    jmin = 3 / 2.;
  }

  for (Spin j = jmin; j <= __jmax; j++) {
    for (Spin m = j * (-1); m <= j; m++) {
      for (Spin n = j * (-1); n <= j; n++) {
        double djmn = 0.;
        for (Spin mm = -1; mm <= Spin(1); mm++) {
          for (Spin nn = -1; nn <= Spin(1); nn++) {
            djmn += Clebsch((double)j - 1, (double)(m - mm), 1, (double)mm,
                            (double)j, (double)m) *
                    Clebsch((double)j - 1, (double)(n - nn), 1, (double)nn,
                            (double)j, (double)n) *
                    __d[j - Spin(1)][m - mm][n - nn] * d1[mm][nn];
          }
        }
        __d[j][m][n] = djmn;
      }
    }
  }
}*/
//_____________________________________________________________________________
/** Fills @a D w/ \f$D^{j}_{m,n}(\alpha,\beta,\gamma)\f$. Uses the recursive
 *  Wigner_d, thus these should be valid to much higher spin than the single
 *  spin version.
 */
/*inline void Wigner_D(const Spin &__jmax, double __alpha, double __beta, double __gamma,
              std::map<Spin, std::map<Spin, std::map<Spin, std::complex<double>>>> &__D) {

  std::complex<double> i(0., 1.);
  std::map<Spin, std::map<Spin, std::map<Spin, double>>> d;
  Wigner_d(__jmax, __beta, d);
  Spin jmin;
  if (d.find(0) != d.end())
    jmin = 0;
  else
    jmin = 1 / 2.;
  for (Spin j = jmin; j <= __jmax; j++) {
    for (Spin m = j * (-1); m <= j; m++) {
      for (Spin n = j * (-1); n <= j; n++)
        __D[j][m][n] =
            exp(-i * ((double)m * __alpha + (double)n * __gamma)) * d[j][m][n];
    }
  }
}*/

} // namespace QFT
} // namespace ComPWA
#endif
