//
//  PhysConstTest.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 24/02/2017.
//
//
#define BOOST_TEST_MODULE                                                      \
  Core /* this can only be define once within the same library ?! */

#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/PhysConst.hpp"
#include "Core/Logging.hpp"

BOOST_AUTO_TEST_SUITE(Core)

BOOST_AUTO_TEST_CASE(XMLInput) {
  ComPWA::Logging log("", boost::log::trivial::severity_level::trace);

  std::stringstream XMLIn;
  XMLIn << "<ParticleList> \
  <Particle Name='gamma'> \
  <Id>22</Id> \
  <Mass Name='mass_gamma'> \
  <Value>0.</Value> \
  </Mass> \
  <Charge>0</Charge> \
  <Spin>1.0</Spin> \
  <Parity>-1</Parity> \
  <Cparity>-1</Cparity> \
  <Gparity>-1</Gparity> \
  <IsoSpin>0</IsoSpin> \
  <IsoSpinZ>0</IsoSpinZ> \
  </Particle> \
  <Particle Name='f0_980'> \
  <Id>9010221</Id> \
  <Mass Name='mass_f0_980'> \
  <Value>0.99</Value> \
  <Fix>1</Fix> \
  <Min>0.5</Min> \
  <Max>1.5</Max> \
  <Error>0</Error> \
  </Mass> \
  <Charge>0</Charge> \
  <Spin>0.0</Spin> \
  <Parity>1</Parity> \
  <Cparity>1</Cparity> \
  <Gparity>1</Gparity> \
  <IsoSpin>0</IsoSpin> \
  <IsoSpinZ>0</IsoSpinZ> \
  <DecayInfo Type='relativisticBreitWigner'> \
  <Width Name='width_f0_980'> \
  <Value>0.05</Value> \
  <Fix>1</Fix> \
  <Min>0.025</Min> \
  <Max>.5</Max> \
  <Error>0</Error> \
  </Width> \
  <MesonRadius Name='radius_f0_980'> \
  <Value>1.5</Value> \
  <Fix>1</Fix> \
  <Min>1.0</Min> \
  <Max>2.0</Max> \
  <Error>0</Error> \
  </MesonRadius> \
  </DecayInfo> \
  </Particle> \
  </ParticleList>";

  boost::property_tree::ptree tr;
  boost::property_tree::xml_parser::read_xml(XMLIn, tr);

  auto inst = ComPWA::PhysConst::CreateInstance(tr);

  auto part = inst->FindParticle("gamma");
  BOOST_CHECK_EQUAL(part.GetMass(), 0.);
  BOOST_CHECK_EQUAL((double)part.GetSpin(), 1.);
  BOOST_CHECK_EQUAL(part.GetParity(), -1);
  BOOST_CHECK_EQUAL(part.GetCparity(), -1);
  BOOST_CHECK_EQUAL(part.GetDecayType(), "stable");

  part = inst->FindParticle(9010221);
  BOOST_CHECK_EQUAL(part.GetMass(), 0.99);
  BOOST_CHECK_EQUAL((double)part.GetSpin(), 0.);
  BOOST_CHECK_EQUAL(part.GetParity(), 1);
  BOOST_CHECK_EQUAL(part.GetCparity(), 1);
  BOOST_CHECK_EQUAL(part.GetDecayType(), "relativisticBreitWigner");
};

BOOST_AUTO_TEST_SUITE_END()
