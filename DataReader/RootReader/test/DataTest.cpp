/*-------------------------------------------------------------------------------
* Copyright (c) 2013 michel.
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the GNU Public License v3.0
* which accompanies this distribution, and is available at
* http://www.gnu.org/licenses/gpl.html
*
* Contributors:
*     michel - initial API and implementation
*-------------------------------------------------------------------------------*/

#include <memory>
#include <vector>
#include <string>

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Number
#include <boost/test/unit_test.hpp>

#include "DataReader/RootReader/RootReader.hpp"

namespace ComPWA {
namespace DataReader {

BOOST_AUTO_TEST_SUITE(RootReaderSuite);

BOOST_AUTO_TEST_CASE(ReadingCheck)
{
  std::string file = "../DataTest-input.root";
  RootReader myReader(file);
  BOOST_CHECK_EQUAL(myReader.GetNEvents(), (unsigned int)100000);

  Event event;
  BOOST_CHECK( event.GetName() != myReader.GetEvent(5).GetName());
}

BOOST_AUTO_TEST_SUITE_END();

} /* namespace DataReader */
} /* namespace ComPWA */

