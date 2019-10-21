// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Define Boost test module
#define BOOST_TEST_MODULE Physics

#include "Physics/PhspVolume.hpp"
#include "Core/Logging.hpp"
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>

using namespace ComPWA::Physics;

// Define Boost test suite
BOOST_AUTO_TEST_SUITE(Physics)

/// Test for KallenFunction. Note that we are not testing the [original Kalen
/// function](https://en.wikipedia.org/wiki/K%C3%A4ll%C3%A9n_function), but the
/// one uses square values as arguments and can be nicely factorised (see
/// [Heron's formula](https://en.wikipedia.org/wiki/Heron%27s_formula)).
BOOST_AUTO_TEST_CASE(KallenFunctionTest) {
  double RelativeTolerance(1e-6);
  // Test values
  BOOST_CHECK_CLOSE(KallenFunction(1., 2., 3.), -8., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(2., 3., 4.), -23., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(3., 4., 5.), -44., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(4., 5., 6.), -71., RelativeTolerance);
  // Test symmetry
  BOOST_CHECK_CLOSE(KallenFunction(1., 2., 3.), -8., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(2., 3., 1.), -8., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(3., 1., 2.), -8., RelativeTolerance);
  BOOST_CHECK_CLOSE(KallenFunction(3., 2., 1.), -8., RelativeTolerance);
}

/// Test application for the calculation of the volume of a phasespace.
/// According to [these lecture
/// notes](http://theory.gsi.de/~knoll/Lecture-notes/1-kinematic.pdf), is should
/// be possible to compute the volume of the phasespace just from the initial CM
/// energy and the masses of the final state particles. For two-particle decays,
/// this can be done analytically, but when there are more particles in the
/// final state, one has to perform a numerical integration.
///
/// Here, we use the decay D0->KsK-K+ as a test case for the three. Test values
/// have been computed with Mathematica.
BOOST_AUTO_TEST_CASE(PhspVolumeTest) {
  // * Define initial state and final state masses
  double s = 1.86483 * 1.86483; // D0
  double m1 = 0.497611;         // K0S
  double m2 = 0.493677;         // K-
  double m3 = 0.493677;         // K+
  // double should_be_vol = 0.541493;
  double should_be_vol = 3.0844;

  // * Test s_lower and s_upper
  std::vector<double> m1m2 = {m1, m2};
  std::vector<double> m1m2m3 = {m1, m2, m3};
  auto range = SRange(s, m1m2m3);
  BOOST_CHECK_CLOSE(range.first, 0.98265, 1e-3);
  BOOST_CHECK_CLOSE(range.second, 1.88006, 1e-3);

  // * Test sample vector generation
  std::vector<double> wantSample = {.7, .9, 1.1, 1.3, 1.5};
  auto nsteps = wantSample.size();
  auto stepSize = wantSample[1] - wantSample[0];
  double lBound = wantSample.front();
  double rBound = wantSample.back();
  auto resultSample = IntegrationSample(lBound, rBound, nsteps);
  double val = lBound;
  for (size_t i = 0; i < resultSample.size(); ++i) {
    BOOST_CHECK_CLOSE(resultSample[i], wantSample[i], 1e-6);
    val += stepSize;
  }
  BOOST_CHECK_CLOSE(resultSample.BinSize, stepSize, 1e-6);

  // * Test PhspVolume for 2 particles
  BOOST_CHECK_CLOSE(PhspVolume(s, m1m2), 5.32194, 1e-4);

  // * Test PhspVolume for 3 particles
  BOOST_CHECK_CLOSE(PhspVolume(s, m1m2m3, 10), should_be_vol, 5);
  BOOST_CHECK_CLOSE(PhspVolume(s, m1m2m3, 100), should_be_vol, 1);
  BOOST_CHECK_CLOSE(PhspVolume(s, m1m2m3, 1000), should_be_vol, 1e-1);
  BOOST_CHECK_CLOSE(PhspVolume(s, m1m2m3, 10000), should_be_vol, 1e-2);
  BOOST_CHECK_CLOSE(PhspVolume(s, m1m2m3, 100000), should_be_vol, 1e-3);

  // * Test PhspVolume for e+e- --> n gamma at sqrt(s) = 100 GeV
  s = 100. * 100.;
  std::vector<double> masses{0., 0.};
  std::cout << PhspVolume(s, masses, 1000) << std::endl;
  masses.push_back(0.);
  std::cout << PhspVolume(s, masses, 1000) << std::endl;
  masses.push_back(0.);
  std::cout << PhspVolume(s, masses, 1000) << std::endl;
  masses.push_back(0.);
  std::cout << PhspVolume(s, masses, 1000) << std::endl;
  // masses.push_back(0.);
  // std::cout << PhspVolume(s, masses, 1000) << std::endl;
  // masses.push_back(0.);
  // std::cout << PhspVolume(s, masses, 1000) << std::endl;
  // masses.push_back(0.);
  // std::cout << PhspVolume(s, masses, 1000) << std::endl;
  // masses.push_back(0.);
  // std::cout << PhspVolume(s, masses, 1000) << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
