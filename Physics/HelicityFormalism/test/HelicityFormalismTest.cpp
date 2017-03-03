//
//  PhysConstTest.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 24/02/2017.
//
//
#define BOOST_TEST_MODULE                                                      \
  HelicityFormalism/* this can only be define once within the same library ?! */

#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/PhysConst.hpp"
#include "Core/Logging.hpp"
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"

BOOST_AUTO_TEST_SUITE(HelicityFormalism)

BOOST_AUTO_TEST_CASE(IncoherentConstruction) {
  ComPWA::Logging log("", boost::log::trivial::severity_level::trace);

  boost::property_tree::ptree tr;
  std::ifstream in("HelicityFormalismTest-input.xml");
  boost::property_tree::xml_parser::read_xml(in, tr);
  
  ComPWA::PhysConst::CreateInstance(tr);

  auto intens = ComPWA::Physics::HelicityFormalism::IncoherentIntensity::Factory(tr);
  
//  auto part = physConst->FindParticle("gamma");
//  BOOST_CHECK_EQUAL(part.GetMass(), 0.);
//  BOOST_CHECK_EQUAL((double)part.GetSpin(), 1.);
//  BOOST_CHECK_EQUAL(part.GetParity(), -1);
//  BOOST_CHECK_EQUAL(part.GetCparity(), -1);
//  BOOST_CHECK_EQUAL(part.GetDecayType(), "stable");
//
//  part = physConst->FindParticle(9010221);
//  BOOST_CHECK_EQUAL(part.GetMass(), 0.99);
//  BOOST_CHECK_EQUAL((double)part.GetSpin(), 0.);
//  BOOST_CHECK_EQUAL(part.GetParity(), 1);
//  BOOST_CHECK_EQUAL(part.GetCparity(), 1);
//  BOOST_CHECK_EQUAL(part.GetDecayType(), "relativisticBreitWigner");
};

BOOST_AUTO_TEST_SUITE_END()
