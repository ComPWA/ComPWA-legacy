// Copyright (c) 2019 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Core

#include "Core/FourMomentum.hpp"
#include "Core/Exceptions.hpp"

#include <boost/test/unit_test.hpp>

#include <memory>
#include <vector>

namespace ComPWA {

BOOST_AUTO_TEST_SUITE(FourMomentumTest)

BOOST_AUTO_TEST_CASE(FourMomentum) {
  ComPWA::FourMomentum p4(1, 2, 3, 4);
  ComPWA::FourMomentum p4B(1, 2, 3, 4);
  BOOST_CHECK_EQUAL(p4, p4B);
  BOOST_CHECK_EQUAL(p4.invariantMassSquared(), 2.0);
  BOOST_CHECK_EQUAL(p4.threeMomentumSquared(), 14.0);

  auto pTot =
      p4 + ComPWA::FourMomentum(std::array<double, 4>{1, 2, 3, 5}); //(2,4,6,9)
  BOOST_CHECK_EQUAL(pTot.invariantMassSquared(), 25.0);

  ComPWA::FourMomentum pSum;
  pSum += p4;
  pSum += p4B;
  BOOST_CHECK_EQUAL(pSum.invariantMassSquared(), 8.0);
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace ComPWA