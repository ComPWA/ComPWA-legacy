// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Data_RootDataIOTest

#include "Data/Root/RootDataIO.hpp"
#include "Core/Logging.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"

#include <boost/test/unit_test.hpp>

#include <memory>
#include <string>

namespace ComPWA {
namespace Data {

BOOST_AUTO_TEST_SUITE(RootData);

BOOST_AUTO_TEST_CASE(SimpleWriteReadCheck) {
  ComPWA::Logging log("", "trace");

  // Generate phsp sample
  std::vector<double> FSMasses = {0.5, 0.5, 0.5};
  ComPWA::Data::Root::RootGenerator gen(1.864, FSMasses);
  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(305896);

  auto sample(ComPWA::Data::generatePhsp(200, gen, RandomGenerator));

  ComPWA::Data::Root::RootDataIO RootIO("trtr");
  RootIO.writeData(sample, "RootReaderTest-output.root");

  auto sampleIn(RootIO.readData("RootReaderTest-output.root"));
  BOOST_CHECK_EQUAL(sample.size(), sampleIn.size());

  std::remove("RootReaderTest-output.root"); // delete file
}

BOOST_AUTO_TEST_SUITE_END();

} // namespace Data
} // namespace ComPWA
