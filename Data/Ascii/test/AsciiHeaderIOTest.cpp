// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Data_AsciiHeaderIOTest

#include "Data/Ascii/AsciiHeaderIO.hpp"
#include "Core/Exceptions.hpp"
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <sstream>
#include <string>

#include <iostream>

namespace ComPWA {
namespace Data {
namespace Ascii {

BOOST_AUTO_TEST_SUITE(AsciiHeaderIOTest);

BOOST_AUTO_TEST_CASE(TestExtractHeaderSuccess) {
  std::string HeaderContent{R"(This
is header
content)"};
  std::stringstream StringStream("<header>" + HeaderContent + "</header>");
  BOOST_CHECK_EQUAL(AsciiHeader::extractHeaderContent(StringStream),
                    HeaderContent + "\n");
}

BOOST_AUTO_TEST_CASE(TestExtractHeaderSuccessOneLine) {
  std::string HeaderContent{"This is header content"};
  std::stringstream StringStream("<header>" + HeaderContent + "</header>");
  BOOST_CHECK_EQUAL(AsciiHeader::extractHeaderContent(StringStream),
                    HeaderContent);
}

BOOST_AUTO_TEST_CASE(TestMissing) {
  std::stringstream MissingEndTag(R"(
    <header>
      PIDs: [1, -5, 23]
      Unit: GeV
      Order: E px py pz
    
    0.97
    -0.00357645 0.0962561   0.0181079   0.170545)");
  std::stringstream FaultyEndTag(R"(
    <header>
      PIDs: [1, -5, 23]
    <header>)");
  std::stringstream FaultyBeginTag(R"(
    <head>
      Unit: GeV
      Order: E px py pz
    </header>)");
  BOOST_CHECK_THROW(AsciiHeader::extractHeaderContent(MissingEndTag),
                    ComPWA::CorruptFile);
  BOOST_CHECK_THROW(AsciiHeader::extractHeaderContent(FaultyEndTag),
                    ComPWA::CorruptFile);
  BOOST_CHECK_THROW(AsciiHeader::extractHeaderContent(FaultyBeginTag),
                    ComPWA::CorruptFile);
  BOOST_CHECK_EQUAL(MissingEndTag.tellg(), 0);
  BOOST_CHECK_EQUAL(FaultyEndTag.tellg(), 0);
  BOOST_CHECK_EQUAL(FaultyBeginTag.tellg(), 0);
}

BOOST_AUTO_TEST_CASE(TestReadYAML) {
  std::string FileName = "Data_AsciiHeaderIOTest-YAMLcorrect.dat";
  std::ifstream File(FileName);
  AsciiHeader Header;
  Header.importYAML(AsciiHeader::extractHeaderContent(File));
  // Test case
  int RemainingWords = 0;
  std::string Word;
  while (File >> Word)
    ++RemainingWords;
  BOOST_CHECK_EQUAL(RemainingWords, 12);
  // Test case
  const std::vector<int> PIDs{-12, 345, 67, -89};
  const std::string Unit = "keV";
  const bool EnergyFirst = true;
  // The test
  BOOST_CHECK_EQUAL(Header.getFinalStatePIDs().size(), PIDs.size());
  for (size_t i = 0; i < PIDs.size(); ++i)
    BOOST_CHECK_EQUAL(Header.getFinalStatePIDs().at(i), PIDs.at(i));
  BOOST_CHECK_EQUAL(Header.getUnit(), Unit);
  BOOST_CHECK_EQUAL(Header.isEnergyFirst(), EnergyFirst);
}

BOOST_AUTO_TEST_CASE(TestWriteReadYAML) {
  // Test case
  const std::vector<int> PIDs{211, 421, -411};
  const std::string Unit = "MeV";
  const bool EnergyFirst = true;
  // Write to file
  std::string FileName = "Data_AsciiHeaderIOTest-YAML.dat";
  AsciiHeader HeaderOut(PIDs, Unit, EnergyFirst);
  std::ofstream FileOut(FileName);
  HeaderOut.dumpToYAML(FileOut);
  FileOut.close();
  // Read from file
  AsciiHeader HeaderIn;
  std::ifstream FileIn(FileName);
  HeaderIn.importYAML(AsciiHeader::extractHeaderContent(FileIn));
  FileIn.close();
  std::remove(FileName.c_str());
  // Test if still the same
  BOOST_CHECK_EQUAL(HeaderIn.getFinalStatePIDs().size(), PIDs.size());
  for (size_t i = 0; i < PIDs.size(); ++i)
    BOOST_CHECK_EQUAL(HeaderIn.getFinalStatePIDs().at(i), PIDs.at(i));
  BOOST_CHECK_EQUAL(HeaderIn.getUnit(), Unit);
  BOOST_CHECK_EQUAL(HeaderIn.isEnergyFirst(), EnergyFirst);
}

BOOST_AUTO_TEST_SUITE_END();

} // namespace Ascii
} // namespace Data
} // namespace ComPWA
