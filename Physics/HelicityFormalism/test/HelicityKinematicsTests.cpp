// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Define Boost test module
#define BOOST_TEST_MODULE HelicityFormalism

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

const std::string HelicityTestParticles = R"####(
<ParticleList>
  <Particle Name='pi0'>
    <Pid>111</Pid>
	<QuantumNumber Class='Spin' Type='Spin' Value='0'/>
	<QuantumNumber Class='Int' Type='Parity' Value='-1'/>
  </Particle>
  <Particle Name='gamma'>
    <Pid>22</Pid>
	<QuantumNumber Class='Spin' Type='Spin' Value='1'/>
	<QuantumNumber Class='Int' Type='Parity' Value='-1'/>
  </Particle>
  <Particle Name='jpsi'>
    <Pid>443</Pid>
	<QuantumNumber Class='Spin' Type='Spin' Value='1'/>
	<QuantumNumber Class='Int' Type='Parity' Value='-1'/>
  </Particle>
</ParticleList>
)####";

// Define Boost test suite (no idea what's the difference to TEST_MODULE)
BOOST_AUTO_TEST_SUITE(HelicityKinematics)

/*! Test application for the calculation of the helicity angle.
 * As example the decay D0->KsK-K+ is used.
 * Some test may be specific to this decay but if all test are passed we can
 * be sure that helicity angles are calculated correctly.
 *
 * A note on numerical precision:
 * We compare all valued using float precision only since otherwise variables
 * calculated using different methods are not completely equal on double
 * precision. Every not an then even with float precision the test fails.
 */
BOOST_AUTO_TEST_CASE(CreateAllSubsystems) {
  ComPWA::Logging log("", "debug");

  // Construct HelicityKinematics from XML tree
  boost::property_tree::ptree tr;
  std::stringstream modelStream;
  // Construct particle list from XML tree
  modelStream << HelicityTestParticles;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto partL = std::make_shared<ComPWA::PartList>();
  ComPWA::ReadParticles(partL, tr);

  // Construct HelicityKinematics by hand
  std::vector<int> finalState, initialState;
  initialState.push_back(443);
  finalState.push_back(22);
  finalState.push_back(111);
  finalState.push_back(111);
  auto kin =
      std::make_shared<ComPWA::Physics::HelicityFormalism::HelicityKinematics>(
          partL, initialState, finalState);
  kin->createAllSubsystems();
  BOOST_CHECK_EQUAL(kin->subSystems().size(), 6);

  finalState.push_back(111);
  auto kin2 =
      std::make_shared<ComPWA::Physics::HelicityFormalism::HelicityKinematics>(
          partL, initialState, finalState);
  kin2->createAllSubsystems();
  BOOST_CHECK_EQUAL(kin2->subSystems().size(), 37);
  finalState.push_back(111);
  auto kin3 =
      std::make_shared<ComPWA::Physics::HelicityFormalism::HelicityKinematics>(
          partL, initialState, finalState);
  kin3->createAllSubsystems();
  BOOST_CHECK_EQUAL(kin3->subSystems().size(), 270);
}

BOOST_AUTO_TEST_SUITE_END()
