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
#ifndef _Utils_H
#define _Utils_H
//_____________________________________________________________________________
/** @file Utils.h
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
#include "Spin.h"
#include "../../include/tensor.h"

using namespace std;
//_____________________________________________________________________________
/// A L,S combination
typedef struct {
  int L;
  Spin S;
} LS;
//_____________________________________________________________________________
/// Returns all valid @a LS combos coupling 1,2 to \f$j^p\f$
vector<LS> GetValidLS(const Spin &__j,int __parity,const Spin &__s1,int __p1,
		      const Spin &__s2,int __p2);
//_____________________________________________________________________________
/**  Calculates \f$ \Gamma (z) \f$
 *   Uses algorithm from C.Lanczos, SIAM Journal of Numerical Analysis 
 *   B1 (1964), 86. 
 *
 *   If \f$ z > 0\f$, returns \f$ \Gamma (z) \f$. If \f$ z = 0\f$ returns 0.
 *   If \f$ z < 0\f$ returns \f$ \frac{-\pi}{|z| \Gamma (|z|) sin(\pi |z|)}\f$.
 */
double Gamma(double __z);
//_____________________________________________________________________________
/** @brief  Calcultates Clebsch-Gordon coeficents.
 *
 * <b> Returns </b>
 *
 *  \f$\left(j_1,m_1;j_2,m_2|J,M\right) \f$ 
 *
 * Note: This function was copied from one written by D.P. Weygand.
 */
double Clebsch(const Spin &__j1,const Spin &__m1,const Spin &__j2,
	       const Spin &__m2,const Spin &__J,const Spin &__M);
//_____________________________________________________________________________
/** Returns \f$d^{j}_{m,n}(\beta)\f$.
 * Calculates the Wigner d-functions. 
 *
 *  This function was copied from one written by D.P. Weygand. It shouldn't be
 *  used for spins higher than about 11/2 or 8.
 */
double Wigner_d(const Spin &__j,const Spin &__m,const Spin &__n,double __beta);
//_____________________________________________________________________________
/** Returns the Wigner D-function \f$D^{j}_{m,n}(\alpha,\beta,\gamma)\f$.
 *  This function uses the single spin Wigner_d function, thus it shares its
 *  limitations regarding higher spins.
 */
inline complex<double> Wigner_D(double __alpha,double __beta,double __gamma,
			const Spin &__j,const Spin &__m,const Spin &__n){
  complex<double> i(0.,1.);
  return exp(-i*(__m*__alpha + __n*__gamma))*Wigner_d(__j,__m,__n,__beta);
}
//_____________________________________________________________________________
/** Fills @a d w/ \f$d^{j}_{m,n}(\beta)\f$ for all spins from 0 to @a jmax 
 *  (integer spin) or 1/2 to @a jmax (half-integer spin). Calculations are 
 *  performed recursively, thus these are valid to much higher spin than the
 *  single-spin version.
 */
void Wigner_d(const Spin &__jmax,double __beta,
	      map<Spin,map<Spin,map<Spin,double> > > &__d);
//_____________________________________________________________________________
/** Fills @a D w/ \f$D^{j}_{m,n}(\alpha,\beta,\gamma)\f$. Uses the recursive
 *  Wigner_d, thus these should be valid to much higher spin than the single
 *  spin version.
 */
void Wigner_D(const Spin &__jmax,double __alpha,double __beta,double __gamma,
	      map<Spin,map<Spin,map<Spin,complex<double> > > > &__D);
//_____________________________________________________________________________
/// Returns \f$ \frac{m\Gamma}{p^2 - m^2 + i m \Gamma} \f$
inline complex<double> BreitWigner(const Vector4<double> &__p4,double __mass,
				   double __width){
  complex<double> i(0.,1.);
  return __mass*__width/(__p4*__p4 - __mass*__mass + i*__mass*__width);
}
//_____________________________________________________________________________
/// Returns \f$ \frac{\Lambda^2 - m^2}{\Lambda^2 - p^2} \f$
inline double MonopoleFormFactor(double __lambda,double __mass,double __p2){
  return (__lambda*__lambda - __mass*__mass)/(__lambda*__lambda - __p2);
}
//_____________________________________________________________________________
/** @brief Calculate the Regge trajectory propagator.
 *
 * @param t Mandelstam variable t
 * @param s Mandelstam variable s
 * @param a Slope of the Regge trajectory
 * @param b Intercept of the Regge trajectory
 * @param spin Spin of the lowest mass member of the family
 * @param sig Signature of the trajectory
 * @param exp_fact Factor to multiple exponential phase term by (defaults to 1)
 *
 * Note: For u-channel, the first argument should be u instead of t.
 *
 * <b> Returns </b> <br>
 *
 * \f$ s^{\alpha(t) - J}\frac{\pi a}{2} \frac{S + F e^{-i\pi\alpha(t)}}
 *     {sin(\pi(\alpha + 1 - J))\Gamma(\alpha + 1 - J)} \f$
 *
 * where \f$\alpha(x) = a\cdot x + b \f$, \f$S\f$ is the signature, \f$F\f$ is
 * the exponential factor and \f$J\f$ is the spin.
 */
complex<double> ReggePropagator(double __t,double __s,double __a,double __b,
				const Spin &__spin,int __sig,
				int __exp_factor = 1);
//_____________________________________________________________________________

/// Gets spin from a string (ex. "1/2" returns a spin of 0.5)
Spin GetSpin(const string &__spin);
//_____________________________________________________________________________

#endif /* _Utils_H */
