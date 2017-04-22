

//
//  PhysConstTest.cpp
//  COMPWA
//
//  Created by Peter Weidenkaff on 24/02/2017.
//
//
#define BOOST_TEST_MODULE                                                      \
  HelicityFormalism /* this can only be define once within the same library ?! \
  */
#include <vector>

#include <locale>
#include <boost/test/unit_test.hpp>
#include <boost/locale/utf.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
//#include <boost/property_tree/xml_parser_writer_settings.hpp>

#include "Core/PhysConst.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Logging.hpp"
#include "Core/RunManager.hpp"
#include "Tools/RootGenerator.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

using namespace ComPWA::Physics::HelicityFormalism;

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <iostream>

BOOST_AUTO_TEST_SUITE(HelicityFormalism)

BOOST_AUTO_TEST_CASE(IncoherentConstruction) {
  
  ComPWA::Logging log("", boost::log::trivial::severity_level::trace);

  boost::property_tree::ptree tr;
  boost::property_tree::xml_parser::read_xml("HelicityFormalismTest-input.xml",
                                             tr);

  ComPWA::PhysConst::CreateInstance(tr);

  //  std::vector<int> finalState, initialState;
  //  initialState.push_back(443);
  //  finalState.push_back(210);
  //  finalState.push_back(210);
  //  finalState.push_back(22);
  //  HelicityKinematics::CreateInstance(initialState, finalState);

  // Create HelicityInstance here
  HelicityKinematics::CreateInstance(tr);
  HelicityKinematics::Instance()->GetPhspVolume();

  // Create amplitude
  auto intens = IncoherentIntensity::Factory(tr);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(123));
  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::RootReader());

  ComPWA::RunManager r;
  r.setGenerator(gen);
  r.setPhspSample(sample);
  r.generatePhsp(2000);

  LOG(info) << "Loop over phsp events....";
  for (auto i : sample->getEvents()) {
    ComPWA::dataPoint p;
    try {
      p = ComPWA::dataPoint(i);
    } catch (std::exception &ex) {
      // Test if events outside the phase space boundaries are generated
      LOG(error) << "Event outside phase space. This should not happen since "
                    "we use a Monte-Carlo sample!";
      BOOST_TEST(false);
      continue;
    }

    double w = intens->Intensity(p);
    i.setWeight(w);
    //    LOG(info) << "point = " << p << " intensity = " << w;
  }

  ComPWA::ParameterList sampleList(sample->getListOfData());
  // Testing function tree
  auto tree = intens->GetTree(sampleList, sampleList, sampleList);
  tree->recalculate();

  std::stringstream printTree;
  printTree << tree;
  LOG(info) << std::endl << printTree.str();

  auto ptout = IncoherentIntensity::Save(intens);
//  if (ptout != tr) {
//    BOOST_CHECK(false);
//    LOG(error) << "Read-in tree and write-out tree are not the same. This is "
//                  "most likely due to an encoding problem but could also "
//                  "point to a bug in reading and writing amplitudes.";
//  }
  
  // Write the property tree to the XML file. Add a line break at the end of each
  // line.
  boost::property_tree::xml_parser::write_xml(
      "HelicityFormalismTest-output.xml", ptout, std::locale(),
      boost::property_tree::xml_writer_make_settings<std::string>(
          ' ', 4, "utf-8"));
};

BOOST_AUTO_TEST_SUITE_END()
