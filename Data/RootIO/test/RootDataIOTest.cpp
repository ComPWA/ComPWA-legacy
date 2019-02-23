// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Data_RootDataIOTest

#include <memory>
#include <string>

#include <boost/test/unit_test.hpp>

#include "Data/DataSet.hpp"
#include "Data/RootIO/RootDataIO.hpp"
#include "Tools/Generate.hpp"
#include "Tools/RootGenerator.hpp"

namespace ComPWA {
namespace Data {

BOOST_AUTO_TEST_SUITE(RootDataIOSuite);

BOOST_AUTO_TEST_CASE(SimpleWriteReadCheck) {
  ComPWA::Logging log("", "trace");

  // Generate phsp sample
  std::vector<double> FSMasses = {0.5, 0.5, 0.5};
  std::shared_ptr<ComPWA::Generator> gen(
      new ComPWA::Tools::RootGenerator(1.864, FSMasses, 305896));

  std::shared_ptr<DataSet> sample(ComPWA::Tools::generatePhsp(200, gen));

  RootDataIO RootIO("trtr");
  RootIO.writeData(sample, "RootReaderTest-output.root");

  std::shared_ptr<DataSet> sampleIn(
      RootIO.readData("RootReaderTest-output.root"));
  BOOST_CHECK_EQUAL(sample->getEventList().size(),
                    sampleIn->getEventList().size());

  std::remove("RootReaderTest-output.root"); // delete file
}

BOOST_AUTO_TEST_SUITE_END();

} // namespace Data
} // namespace ComPWA
