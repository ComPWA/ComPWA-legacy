// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

/// \file
/// Function to calculate the available phase space volume of arbitrary decays.
/// The algorithm follows the recursive integral computation outlined in [this
/// paper](http://theory.gsi.de/~knoll/Lecture-notes/1-kinematic.pdf), see
/// pp.&nbsp;6–7.

#ifndef PhspVolume_h
#define PhspVolume_h

#include <utility>
#include <vector>

namespace ComPWA {
namespace Physics {

/// Compute phasespace volume of momentum space for an arbitrary number of
/// particles in the final state using Riemann integration.
/// @return A pair: first value is the volume, second is the error (currently
/// set to `0.`)
/// @todo Implement errors (second member of the pair).
/// @todo Algorithm might be improved with [Simpson's
/// rule](https://en.wikipedia.org/wiki/Simpson%27s_rule), because we integrate
/// over a function that is polynomial in the limit \f$m_i\rightarrow 0\f$
std::pair<double, double> PhspVolume(double s, std::vector<double> &FSMasses,
                                     std::size_t SampleSize = 1000);

/// Original [Källén
/// function](https://en.wikipedia.org/wiki/K%C3%A4ll%C3%A9n_function), that is,
/// not having square values in its argument. We use this function instead of
/// the one that can be factorised (see [Heron's
/// formula](https://en.wikipedia.org/wiki/Heron%27s_formula)), because we need
/// to enter \f$s\f$ without taking its square root.
double KallenFunction(double x, double y, double z);

/// Phase space element for a two particle decay. An analytic solution exists
/// only for the volume of the phasespace of two-particle decays.
double PhspVolumeTwoParticles(double s, double m1, double m2);

} // namespace Physics
} // namespace ComPWA

#endif
