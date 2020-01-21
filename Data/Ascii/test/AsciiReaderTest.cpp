// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Data_AsciiReaderTest

#include "Data/Ascii/AsciiReader.hpp"
#include "Core/Logging.hpp"
#include "Data/EvtGen/EvtGenGenerator.hpp"
#include "Data/Generate.hpp"

#include <boost/test/unit_test.hpp>
#include <string>

#include <fstream>

namespace ComPWA {
namespace Data {
namespace Ascii {

std::vector<ComPWA::Event>
generateSample(std::size_t NumberOfEvents,
               const ComPWA::FourMomentum &CMSP4 = 1.85,
               const std::vector<double> &masses = {.5, .5, .5}) {
  auto EventGenerator = ComPWA::Data::EvtGen::EvtGenGenerator(CMSP4, masses);
  ComPWA::StdUniformRealGenerator RandomGenerator(1234);
  std::vector<ComPWA::Event> Events(NumberOfEvents);
  for (auto &Event : Events) {
    Event = EventGenerator.generate(RandomGenerator);
    Event.Weight = RandomGenerator();
  }
  return Events;
}

BOOST_AUTO_TEST_SUITE(AsciiData);

BOOST_AUTO_TEST_CASE(TestMomentumEnergyOrder) {
  std::list<std::string> TestFiles{
      "Data_AsciiReaderTest-CorrectWeightedMomE.dat",
      "Data_AsciiReaderTest-CorrectWeightedEmom.dat"};
  for (const auto &FileName : TestFiles) {
    auto Events = readData(FileName, 10);
    BOOST_CHECK_EQUAL(Events.size(), 3);

    const auto &Event = Events.front();
    BOOST_CHECK_EQUAL(Event.Weight, .77);

    auto Particle = Event.ParticleList.at(0);
    BOOST_CHECK_EQUAL(Particle.pid(), 123);
    BOOST_CHECK_EQUAL(Particle.fourMomentum().e(), 5.);
    BOOST_CHECK_EQUAL(Particle.fourMomentum().px(), .543);
    BOOST_CHECK_EQUAL(Particle.fourMomentum().py(), .2345);
    BOOST_CHECK_EQUAL(Particle.fourMomentum().pz(), 1.);

    Particle = Event.ParticleList.at(2);
    BOOST_CHECK_EQUAL(Particle.pid(), -123);
    BOOST_CHECK_EQUAL(Particle.fourMomentum().e(), 9.);
    BOOST_CHECK_EQUAL(Particle.fourMomentum().px(), .85434);
    BOOST_CHECK_EQUAL(Particle.fourMomentum().py(), .564);
    BOOST_CHECK_EQUAL(Particle.fourMomentum().pz(), .923);
  }
}

BOOST_AUTO_TEST_CASE(TestOverwrite) {
  const char *Filename = "Data_AsciiReaderTest-TestOverwrite.dat";
  auto Events1 = generateSample(3);
  writeData(Events1, Filename);
  auto Events2 = generateSample(4);
  writeData(Events2, Filename, true);
  auto ImportedEvents = readData(Filename);
  BOOST_CHECK_EQUAL(ImportedEvents.size(), Events1.size() + Events2.size());
  std::remove(Filename);
}

BOOST_AUTO_TEST_SUITE_END();

} // namespace Ascii
} // namespace Data
} // namespace ComPWA
