// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE Core

#include "Core/Properties.hpp"
#include "Core/Logging.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>
#include <sstream>

BOOST_AUTO_TEST_SUITE(Core)

BOOST_AUTO_TEST_CASE(XMLInput) {
  ComPWA::Logging log("", "trace");

  std::stringstream XMLIn;
  XMLIn << "<ParticleList> \
  <Particle Name='gamma'> \
  <Pid>22</Pid> \
  <Parameter Type='Mass' Name='mass_gamma'> \
  <Value>0.</Value> \
  </Parameter> \
  <QuantumNumber Class='Spin' Type='Spin' Value='1.0'/> \
  <QuantumNumber Class='Int' Type='Charge' Value='0'/> \
  <QuantumNumber Class='Int' Type='Parity' Value='-1'/> \
  <QuantumNumber Class='Int' Type='Cparity' Value='-1'/> \
  <QuantumNumber Class='Int' Type='Gparity' Value='-1'/> \
  </Particle> \
  <Particle Name='f0_980'> \
  <Pid>9010221</Pid> \
  <Parameter Type='Mass' Name='mass_f0_980'> \
  <Value>0.99</Value> \
  <Fix>1</Fix> \
  <Min>0.5</Min> \
  <Max>1.5</Max> \
  <Error>0</Error> \
  </Parameter> \
  <QuantumNumber Class='Spin' Type='Spin' Value='0.0'/> \
  <QuantumNumber Class='Int' Type='Charge' Value='0'/> \
  <QuantumNumber Class='Int' Type='Parity' Value='1'/> \
  <QuantumNumber Class='Int' Type='Cparity' Value='1'/> \
  <QuantumNumber Class='Int' Type='Gparity' Value='1'/> \
  <DecayInfo Type='relativisticBreitWigner'> \
  <Parameter Type='Width' Name='width_f0_980'> \
  <Value>0.05</Value> \
  <Fix>1</Fix> \
  <Min>0.025</Min> \
  <Max>.5</Max> \
  <Error>0</Error> \
  </Parameter> \
  <Parameter Type='MesonRadius' Name='radius_f0_980'> \
  <Value>1.5</Value> \
  <Fix>1</Fix> \
  <Min>1.0</Min> \
  <Max>2.0</Max> \
  <Error>0</Error> \
  </Parameter> \
  </DecayInfo> \
  </Particle> \
  </ParticleList>";

  auto partL = ComPWA::readParticles(XMLIn);

  auto part = ComPWA::findParticle(partL, "gamma");
  BOOST_CHECK_EQUAL(part.getMass().Value, 0.);
  BOOST_CHECK_EQUAL(part.getQuantumNumber<double>("Spin"), 1.);
  BOOST_CHECK_EQUAL(part.getQuantumNumber<int>("Parity"), -1);
  BOOST_CHECK_EQUAL(part.getQuantumNumber<int>("Cparity"), -1);
  BOOST_CHECK_EQUAL(part.getDecayType(), "Stable");

  part = ComPWA::findParticle(partL, 9010221);
  BOOST_CHECK_EQUAL(part.getMass().Value, 0.99);
  BOOST_CHECK_EQUAL(part.getQuantumNumber<double>("Spin"), 0.);
  BOOST_CHECK_EQUAL(part.getQuantumNumber<int>("Parity"), 1);
  BOOST_CHECK_EQUAL(part.getQuantumNumber<int>("Cparity"), 1);
  BOOST_CHECK_EQUAL(part.getDecayType(), "relativisticBreitWigner");
};

BOOST_AUTO_TEST_SUITE_END()