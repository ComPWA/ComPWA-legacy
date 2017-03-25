//
//  PhysConstTest.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 24/02/2017.
//
//
#define BOOST_TEST_MODULE                                                      \
  HelicityFormalism/* this can only be define once within the same library ?! */
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/PhysConst.hpp"
#include "Core/Logging.hpp"
#include "Core/RunManager.hpp"
#include "Tools/RootGenerator.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

BOOST_AUTO_TEST_SUITE(HelicityFormalism)

BOOST_AUTO_TEST_CASE(IncoherentConstruction) {
  ComPWA::Logging log("", boost::log::trivial::severity_level::trace);

  boost::property_tree::ptree tr;
  std::ifstream in("HelicityFormalismTest-input.xml");
  boost::property_tree::xml_parser::read_xml(in, tr);
  
  ComPWA::PhysConst::CreateInstance(tr);
  
  std::vector<int> finalState, initialState;
  initialState.push_back(443);
  finalState.push_back(210);
  finalState.push_back(210);
  finalState.push_back(22);
  
  //Create HelicityInstance here
  ComPWA::Physics::HelicityFormalism::HelicityKinematics::CreateInstance(initialState,finalState);
  
  //Create amplitude
  auto intens = ComPWA::Physics::HelicityFormalism::IncoherentIntensity::Factory(tr);
  
  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen( new ComPWA::Tools::RootGenerator(0) );
  std::shared_ptr<ComPWA::DataReader::Data> sample( new ComPWA::DataReader::RootReader() );
  ComPWA::RunManager r;
  r.setGenerator(gen);
  r.setPhspSample(sample);
  r.generatePhsp(10);
  auto dataP = sample->getDataPoints();
  
  
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
