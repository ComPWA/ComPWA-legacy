// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Data_AsciiDataIOTest

#include "Data/Ascii/AsciiDataIO.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Logging.hpp"
#include "Data/EvtGen/EvtGenGenerator.hpp"
#include "Data/Generate.hpp"

#include <boost/test/unit_test.hpp>
#include <string>

#include <fstream>

namespace ComPWA {
namespace Data {
namespace Ascii {

EventCollection generateSample(std::size_t NumberOfEvents,
                               const std::vector<pid> &Pids,
                               const std::vector<double> &masses,
                               const ComPWA::FourMomentum &CMSP4 = 1.85) {
  auto EventGenerator =
      ComPWA::Data::EvtGen::EvtGenGenerator(CMSP4, masses, Pids);
  ComPWA::StdUniformRealGenerator RandomGenerator(1234);
  return EventGenerator.generate(NumberOfEvents, RandomGenerator);
}

BOOST_AUTO_TEST_SUITE(AsciiData);

BOOST_AUTO_TEST_CASE(TestCorrectUnweighted) {
  std::string FileName = "Data_AsciiDataIOTest-CorrectUnweighted.dat";
  auto EvtList = readData(FileName);
  BOOST_CHECK_EQUAL(EvtList.Events.size(), 3);

  BOOST_CHECK_EQUAL(EvtList.Pids.at(0), 123);
  BOOST_CHECK_EQUAL(EvtList.Pids.at(2), -123);

  const auto &Event = EvtList.Events.front();
  BOOST_CHECK_EQUAL(Event.Weight, 1.);

  auto FourVector = Event.FourMomenta.at(0);
  BOOST_CHECK_EQUAL(FourVector.e(), 5.);
  BOOST_CHECK_EQUAL(FourVector.px(), .543);
  BOOST_CHECK_EQUAL(FourVector.py(), .2345);
  BOOST_CHECK_EQUAL(FourVector.pz(), 1.);

  FourVector = Event.FourMomenta.at(2);
  BOOST_CHECK_EQUAL(FourVector.e(), 9.);
  BOOST_CHECK_EQUAL(FourVector.px(), .85434);
  BOOST_CHECK_EQUAL(FourVector.py(), .564);
  BOOST_CHECK_EQUAL(FourVector.pz(), .923);
}

BOOST_AUTO_TEST_CASE(TestMomentumEnergyOrder) {
  std::list<std::string> TestFiles{
      "Data_AsciiDataIOTest-CorrectWeightedMomE.dat",
      "Data_AsciiDataIOTest-CorrectWeightedEmom.dat"};
  for (const auto &FileName : TestFiles) {
    auto EvtList = readData(FileName, 10);
    BOOST_CHECK_EQUAL(EvtList.Events.size(), 3);

    BOOST_CHECK_EQUAL(EvtList.Pids.at(0), 123);
    BOOST_CHECK_EQUAL(EvtList.Pids.at(2), -123);

    const auto &Event = EvtList.Events.front();
    BOOST_CHECK_EQUAL(Event.Weight, .77);

    auto FourVector = Event.FourMomenta.at(0);
    BOOST_CHECK_EQUAL(FourVector.e(), 5.);
    BOOST_CHECK_EQUAL(FourVector.px(), .543);
    BOOST_CHECK_EQUAL(FourVector.py(), .2345);
    BOOST_CHECK_EQUAL(FourVector.pz(), 1.);

    FourVector = Event.FourMomenta.at(2);
    BOOST_CHECK_EQUAL(FourVector.e(), 9.);
    BOOST_CHECK_EQUAL(FourVector.px(), .85434);
    BOOST_CHECK_EQUAL(FourVector.py(), .564);
    BOOST_CHECK_EQUAL(FourVector.pz(), .923);
  }
}

BOOST_AUTO_TEST_CASE(TestOverwrite) {
  const std::string FileName = "Data_AsciiDataIOTest-TestOverwrite.dat";
  std::vector<pid> Pids = {1, 2, 3};
  auto Sample1 = generateSample(12, Pids, {0.5, 0.5, 0.5});
  writeData(Sample1, FileName);
  auto Sample2 = generateSample(17, Pids, {0.5, 0.5, 0.5});
  writeData(Sample2, FileName, false);
  auto ImportedEvents = readData(FileName);
  BOOST_CHECK_EQUAL(ImportedEvents.Events.size(),
                    Sample1.Events.size() + Sample2.Events.size());
  for (size_t i = 0; i < ImportedEvents.Events.size(); ++i) {
    const auto &EventIn = ImportedEvents.Events[i];
    Event EventOut;
    if (i < Sample1.Events.size())
      EventOut = Sample1.Events[i];
    else
      EventOut = Sample2.Events[i - Sample1.Events.size()];
    BOOST_CHECK_CLOSE_FRACTION(EventIn.Weight, EventOut.Weight, 1e-5);
    assert(EventIn.FourMomenta.size() == EventOut.FourMomenta.size());
    for (size_t j = 0; j < EventIn.FourMomenta.size(); ++j)
      BOOST_CHECK_CLOSE_FRACTION(EventIn.FourMomenta.at(j).invariantMass(),
                                 EventOut.FourMomenta.at(j).invariantMass(),
                                 1e-5);
  }
  std::remove(FileName.c_str());
}

BOOST_AUTO_TEST_SUITE_END();

} // namespace Ascii
} // namespace Data
} // namespace ComPWA
