//
//  PhysConstTest.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 24/02/2017.
//
//
#define BOOST_TEST_MODULE Core /* this can only be define once within the same library ?! */

#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/PhysConst.hpp"

BOOST_AUTO_TEST_SUITE( Core )

BOOST_AUTO_TEST_CASE( xmlInput )
{
  std::stringstream partXML;
  partXML << "<ParticleList> \
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
  </ParticleList>";
  
  boost::property_tree::ptree tr;
  boost::property_tree::xml_parser::read_xml( partXML, tr );
  
  auto inst = ComPWA::PhysConst::CreateInstance(tr);
  
  auto part = inst->FindParticle("gamma");
  BOOST_CHECK_EQUAL( part.GetMass() , 0. );
  BOOST_CHECK_EQUAL( (double)part.GetSpin() , 1. );
  BOOST_CHECK_EQUAL( part.GetParity() , -1 );
  BOOST_CHECK_EQUAL( part.GetCparity() , -1 );
  
  
  
  //  auto inst = ComPWA::PhysConst::createInstance("particles.xml");
  //  const_string cs1( "test_string" );

  //  BOOST_CHECK_THROW( cs1.at( cs1.length() ), std::out_of_range );
};
BOOST_AUTO_TEST_SUITE_END()
