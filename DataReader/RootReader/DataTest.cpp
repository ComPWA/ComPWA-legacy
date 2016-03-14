#-------------------------------------------------------------------------------
# Copyright (c) 2013 michel.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/gpl.html
# 
# Contributors:
#     michel - initial API and implementation
#-------------------------------------------------------------------------------
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Number
#include <boost/test/unit_test.hpp>
#include "DataReader/RootReader/RootReader.hpp"
#include "Core/PWAEvent.hpp"
#include <memory>
#include <vector>
#include <string>

namespace ComPWA {
namespace DataReader {
namespace RootReader {

BOOST_AUTO_TEST_SUITE(RootReaderSuite);

BOOST_AUTO_TEST_CASE(ReadingCheck)
{
  std::string file = "test/2Part-4vecs.root";
  RootReader myReader(file, false);
  BOOST_CHECK_EQUAL(myReader.getNEvents(), (unsigned int)100000);

  Event event;
  BOOST_CHECK(myReader.getEvent(5, event));
}

BOOST_AUTO_TEST_SUITE_END();

} /* namespace RootReader */
} /* namespace DataReader */
} /* namespace ComPWA */

