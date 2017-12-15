// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <memory>
#include <vector>
#include <string>
#include <cstdio>

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Number
#include <boost/test/unit_test.hpp>

#include "DataReader/RootReader/RootReader.hpp"
#include "Tools/Generate.hpp"
#include "Tools/RootGenerator.hpp"

namespace ComPWA {
namespace DataReader {

BOOST_AUTO_TEST_SUITE(RootReaderSuite);

BOOST_AUTO_TEST_CASE(WriteReadCheck) {
  ComPWA::Logging log("", boost::log::trivial::severity_level::trace);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(
      new ComPWA::Tools::RootGenerator(1.864, 0.5, 0.5, 0.5, 305896));

  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::RootReader());

  //ComPWA::Tools::RunManager r;
  //r.SetGenerator(gen);
  //r.SetPhspSample(sample);
  //r.GeneratePhsp(200);
  ComPWA::Tools::GeneratePhsp(200, gen, sample);

  sample->writeData("RootReaderTest-output.root", "trtr");

  std::shared_ptr<ComPWA::DataReader::Data> sampleIn(
      new ComPWA::DataReader::RootReader("RootReaderTest-output.root", "trtr",
                                         -1));
  BOOST_CHECK_EQUAL(sample->numEvents(), sampleIn->numEvents());
  
  std::remove("RootReaderTest-output.root"); // delete file
}

BOOST_AUTO_TEST_SUITE_END();

} /* namespace DataReader */
} /* namespace ComPWA */
