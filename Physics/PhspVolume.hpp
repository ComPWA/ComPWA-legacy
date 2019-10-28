// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Function to calculate the available phase space volume of arbitrary decays. The algorithm follows the recursive integral computation outlined in [this paper](http://theory.gsi.de/~knoll/Lecture-notes/1-kinematic.pdf), see pp.&nbsp;6–7.
///


#ifndef PhspVolume_h
#define PhspVolume_h

#include <cmath>
#include <vector>

namespace ComPWA {
namespace Tools {

/// Kallen function.
/// (https://en.wikipedia.org/wiki/Källén_function)
inline double KallenFunction(double x, double y, double z){
  return (x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * z * x);
}

/// Threshold function. Can be interpreted as the particle velocity at a given
/// center-of-mass energy \p s in case of equal masses (\p m1 = \p m2).
/// Vanishes below threshold.
inline double ThresholdFunction(double s, double m1, double m2){
  return std::sqrt(KallenFunction(s, m1 * m1, m2 * m2)) / s;
}

/// Phase space element for a two particle decay.
inline double PhspTwoParticles(double s, double m1, double m2){
//  return M_PI / 2 * ThresholdFunction(s, m1, m2); //Tord Riemann
  return ThresholdFunction(s, m1, m2) / (8 * M_PI); //Hitoshi Murayama
}

inline double PhspNParaticles(double s, std::vector<double> masses) {
  double vol;
  if (masses.size() == 2)
    vol = PhspTwoParticles(s, masses.at(0), masses.at(1));
  return vol;
}

} // ns::Tools
} // ns::ComPWA

#endif
