// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Core

#include "Core/Spin.hpp"

#include <boost/test/unit_test.hpp>

namespace ComPWA {

BOOST_AUTO_TEST_SUITE(SpinTest);

void checkFractionsEqual(ComPWA::Fraction f1, ComPWA::Fraction f2) {
  BOOST_CHECK_EQUAL(f1.getNumerator(), f2.getNumerator());
  BOOST_CHECK_EQUAL(f1.getDenominator(), f2.getDenominator());
}

BOOST_AUTO_TEST_CASE(FractionConstruction) {
  checkFractionsEqual(ComPWA::Fraction(1, 2), ComPWA::Fraction(0.5));
  checkFractionsEqual(ComPWA::Fraction(-1, 2), ComPWA::Fraction(-0.5));
  checkFractionsEqual(ComPWA::Fraction(3, 2), ComPWA::Fraction(1.5));
  checkFractionsEqual(ComPWA::Fraction(10, 2), ComPWA::Fraction(5.001));
  checkFractionsEqual(ComPWA::Fraction(8, 16), ComPWA::Fraction(0.5));

  BOOST_CHECK_THROW(ComPWA::Fraction(1, 0), std::exception);
  BOOST_CHECK_THROW(ComPWA::Fraction(-1, 0), std::exception);
  BOOST_CHECK_NO_THROW(ComPWA::Fraction(0.0)); // 0.0 is constructed as 0/1
}

BOOST_AUTO_TEST_CASE(FractionSubtraction) {
  checkFractionsEqual(ComPWA::Fraction(1, 2) - ComPWA::Fraction(0.5), 0.0);
  checkFractionsEqual(ComPWA::Fraction(1, 2) - ComPWA::Fraction(-0.5), 1.0);
  checkFractionsEqual(ComPWA::Fraction(1, 1) - ComPWA::Fraction(0.5), 0.5);
  checkFractionsEqual(ComPWA::Fraction(2, 1) - ComPWA::Fraction(-0.5), 2.5);
}

BOOST_AUTO_TEST_CASE(SpinConstruction) {
  BOOST_CHECK_NO_THROW(ComPWA::Spin(1.0, 0.0));
  BOOST_CHECK_THROW(ComPWA::Spin(-1.0, 0.0), ComPWA::BadParameter);
  BOOST_CHECK_THROW(ComPWA::Spin(1.0, 2.0), ComPWA::BadParameter);
  BOOST_CHECK_THROW(ComPWA::Spin(1.0, -1.5), ComPWA::BadParameter);
  BOOST_CHECK_THROW(ComPWA::Spin(0.5, 0.0), ComPWA::BadParameter);
  BOOST_CHECK_NO_THROW(ComPWA::Spin(0.5, 0.5));
  BOOST_CHECK_THROW(ComPWA::Spin(-0.5, 0.5), ComPWA::BadParameter);
}

BOOST_AUTO_TEST_CASE(SpinValue) {}

BOOST_AUTO_TEST_SUITE_END();

} /* namespace ComPWA */
