//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

#ifndef PHYSICS_DYNAMICALDECAYFUNCTIONS_KINEMATICS_HPP_
#define PHYSICS_DYNAMICALDECAYFUNCTIONS_KINEMATICS_HPP_

#include <complex>

namespace DynamicalFunctions {

class Kinematics {
public:
  Kinematics();
  virtual ~Kinematics() {
  }
  ;

  /** Convert width of resonance to coupling
   *
   * Implementation of Eq.47-21 of PDG2014. Only valid for narrow, isolated resonances.
   * @param mSq invariant mass
   * @param mR mass of resonance
   * @param width width of resonance in channel [a,b]
   * @param ma mass of particle a
   * @param mb mass of particle b
   * @return
   */
  static std::complex<double> widthToCoupling(double mSq, double mR,
      double width, double ma, double mb, double spin, double mesonRadius);
  /** Convert coupling to width
   *
   * Convert coupling to channel (@ma,@mb) to partial width. Only valid for narrow, isolated resonances.
   * Implementation of inverted Eq.47-21 of PDG2014.
   * @param mSq invariant mass
   * @param mR mass of resonance
   * Implementation of inverted Eq.47-21 of PDG2014. Only valid for narrow, isolated resonances.
   * @param mSq invariant mass
   * @param mR mass of resonance
   * @param g coupling to channel [a,b]
   * @param ma mass of particle a
   * @param mb mass of particle b
   * @return
   */
  static std::complex<double> couplingToWidth(double mSq, double mR, double g,
      double ma, double mb, double spin, double mesonRadius);
  /** Calculate Break-up momentum squared
   *
   * Calculate Break-up momentum at energy @sqrtS for particles with masses @ma and @mb.
   * From PDG2014 Eq.46-20a. Below threshold the function is analytically continued.
   * @param sqrtS center-of-mass energy
   * @param ma mass particle A
   * @param mb mass particle B
   * @return |break-up momentum|
   */
  static double qSqValue(double sqrtS, double ma, double mb);
  /** Calculate Break-up momentum
   *
   * Calculate Break-up momentum at energy @sqrtS for particles with masses @ma and @mb.
   * From PDG2014 Eq.46-20a. Below threshold the function is analytically continued.
   * @param sqrtS center-of-mass energy
   * @param ma mass particle A
   * @param mb mass particle B
   * @return |break-up momentum|
   */
  static std::complex<double> qValue(double sqrtS, double ma, double mb);
  /** Two body phsp factor
   *
   * From PDG2014 Eqn.47-2
   * @param sqrtS invariant mass of particles A and B
   * @param ma Mass of particle A
   * @param mb Mass of particle B
   * @return
   */
  static std::complex<double> phspFactor(double sqrtS, double ma, double mb);
  static double FormFactor(double sqrtS, double ma, double mb, double spin,
      double mesonRadius);
};

}
#endif /* PHYSICS_DYNAMICALDECAYFUNCTIONS_KINEMATICS_HPP_ */
