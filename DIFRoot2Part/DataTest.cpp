#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Number
#include <boost/test/unit_test.hpp>
#include "DIFRootReader.hpp"
#include "PWAEvent.hpp"
#include <memory>
#include <vector>
#include <string>

BOOST_AUTO_TEST_SUITE(RootReaderSuite);

BOOST_AUTO_TEST_CASE(ReadingCheck)
{
  std::string file = "test/2Part-4vecs.root";
  DIFRootReader myReader(file, false);
  BOOST_CHECK_EQUAL(myReader.getNEvents(), (unsigned int)100000);

  PWAEvent event;
  BOOST_CHECK(myReader.getEvent(5, event));
}

BOOST_AUTO_TEST_SUITE_END();
