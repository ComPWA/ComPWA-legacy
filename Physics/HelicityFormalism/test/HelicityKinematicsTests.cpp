// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// Define Boost test module
#define BOOST_TEST_MODULE HelicityFormalism

#include "Core/Logging.hpp"
#include "Data/DataSet.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

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

BOOST_AUTO_TEST_CASE(CreateAllSubsystems) {
  ComPWA::Logging log("", "debug");

  std::stringstream modelStream;
  modelStream << HelicityTestParticles;
  auto partL = ComPWA::readParticles(modelStream);

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
  auto Dataset = kin->convert({});
  size_t MassVariables = 1 + 3;
  size_t AngleVariables = (3 + 3) * 2;
  BOOST_CHECK_EQUAL(Dataset.Data.size(), MassVariables + AngleVariables);

  finalState.push_back(111);
  auto kin2 =
      std::make_shared<ComPWA::Physics::HelicityFormalism::HelicityKinematics>(
          partL, initialState, finalState);
  kin2->createAllSubsystems();
  Dataset = kin2->convert({});
  MassVariables = 1 + 4 + 6;
  AngleVariables = (12 + 12 + 4) * 2 + (6 + 3) * 2;
  BOOST_CHECK_EQUAL(Dataset.Data.size(), MassVariables + AngleVariables);
  finalState.push_back(111);
  auto kin3 =
      std::make_shared<ComPWA::Physics::HelicityFormalism::HelicityKinematics>(
          partL, initialState, finalState);
  kin3->createAllSubsystems();
  Dataset = kin3->convert({});
  MassVariables = 1 + 5 + 10 + 10;
  // The angles combinations are grouped by topology. The two factors subtracted
  // at the end are undoing double counting between topologies (There is an
  // overlap between the first and third group: the graph parts most bottom or
  // top are identical)
  AngleVariables =
      ((30 + 15 + 5) + (30 + 30 + 10 + 10) + (30 + 60 + 20 + 5) - 30 - 5) * 2;
  BOOST_CHECK_EQUAL(Dataset.Data.size(), MassVariables + AngleVariables);
}

BOOST_AUTO_TEST_SUITE_END()
