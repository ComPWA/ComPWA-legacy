// Copyright (c) 2019 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE GoodnessOfFitTest

#include "Tools/GoodnessOfFit.hpp"
#include "Core/Logging.hpp"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(ToolTests)

BOOST_AUTO_TEST_CASE(ICtest) {
  ComPWA::Logging Log("trace", "");

  BOOST_CHECK_EQUAL(calculateAIC(100, 12, 3), 209);
  BOOST_CHECK_CLOSE(calculateBIC(100, 12, 3), 207.4547199, 0.001);
}

BOOST_AUTO_TEST_SUITE_END()
