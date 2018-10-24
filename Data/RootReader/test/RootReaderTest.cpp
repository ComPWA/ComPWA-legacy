// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Number
#include <boost/test/unit_test.hpp>

#include "Data/RootReader/RootReader.hpp"
#include "Tools/Generate.hpp"
#include "Tools/RootGenerator.hpp"

namespace ComPWA {
namespace Data {

BOOST_AUTO_TEST_SUITE(RootReaderSuite);

BOOST_AUTO_TEST_CASE(WriteReadCheck) {
  ComPWA::Logging log("", "trace");

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(
      new ComPWA::Tools::RootGenerator(1.864, 0.5, 0.5, 0.5, 305896));

  std::shared_ptr<ComPWA::Data::Data> sample(
      ComPWA::Tools::generatePhsp(200, gen));

  RootReader RootIO("trtr");
  RootIO.writeData(sample, "RootReaderTest-output.root");

  std::shared_ptr<ComPWA::Data::Data> sampleIn(
      RootIO.readData("RootReaderTest-output.root"));
  BOOST_CHECK_EQUAL(sample->numEvents(), sampleIn->numEvents());

  std::remove("RootReaderTest-output.root"); // delete file
}

BOOST_AUTO_TEST_SUITE_END();

} // namespace DataReader
} // namespace ComPWA
