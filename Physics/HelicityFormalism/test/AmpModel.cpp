// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// This can only be define once within the same library ?!
#define BOOST_TEST_MODULE HelicityFormalism

#include <vector>
#include <locale>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/locale/utf.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

#include "Core/Properties.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Logging.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/Generate.hpp"
#include "DataReader/Data.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

// Amplitude model
#include "Physics/HelicityFormalism/test/AmpModelTest.hpp"

using namespace ComPWA;
using namespace ComPWA::Physics::HelicityFormalism;

BOOST_AUTO_TEST_SUITE(HelicityFormalism)

BOOST_AUTO_TEST_CASE(KinematicsConstructionFromXML) {
  ComPWA::Logging log("", boost::log::trivial::severity_level::trace);

  boost::property_tree::ptree tr;
  std::stringstream modelStream;
  // Construct particle list from XML tree
  modelStream << HelicityTestParticles;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, tr);

  modelStream.clear();
  tr = boost::property_tree::ptree();
  // Construct Kinematics from XML tree
  modelStream << HelicityTestKinematics;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto kin = std::make_shared<HelicityKinematics>(
      partL, tr.get_child("HelicityKinematics"));

  BOOST_CHECK_EQUAL(kin->phspVolume(), 0.123);
  BOOST_CHECK_EQUAL(kin->finalState().size(), 3);
  BOOST_CHECK_EQUAL(kin->initialState().size(), 1);
}

BOOST_AUTO_TEST_CASE(ConstructionFromXML) {
  // Due to the structure of Boost.UnitTest the instances already exist from
  // previous test
  //  ComPWA::Logging log("", boost::log::trivial::severity_level::trace);

  boost::property_tree::ptree tr;
  std::stringstream modelStream;
  // Construct particle list from XML tree
  modelStream << HelicityTestParticles;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, tr);

  modelStream.clear();
  tr = boost::property_tree::ptree();
  // Construct Kinematics from XML tree
  modelStream << HelicityTestKinematics;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto kin = std::make_shared<HelicityKinematics>(
      partL, tr.get_child("HelicityKinematics"));

  modelStream.clear();
  tr = boost::property_tree::ptree();
  modelStream << HelicityTestModel;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  // Create amplitude from property_tree
  auto intens = ComPWA::Physics::IncoherentIntensity::Factory(
      partL, kin, tr.get_child("IncoherentIntensity"));

  // Save amplitude to property_tree
  boost::property_tree::ptree ptout;
  ptout.add_child("IncoherentIntensity",
                  ComPWA::Physics::IncoherentIntensity::Save(intens));

  //  if (ptout != tr) {
  //    BOOST_CHECK(false);
  //    LOG(error)
  //        << "Read-in tree and write-out tree are not the same. This is"
  //           "most likely due to an encoding problem but could also "
  //           "point to a bug in reading and writing amplitudes. Check input"
  //           "and output files carefully.";
  //  }

  // Write the property tree to the XML file. Add a line break at the end of
  // each line.
  boost::property_tree::xml_parser::write_xml("AmpModel-output.xml", ptout,
                                              std::locale());

  std::remove("AmpModel-output.xml"); // delete file
  // Compile error for some boost/compiler versions
  //  boost::property_tree::xml_parser::write_xml(
  //      "../HelicityFormalismTest-output.xml", ptout, std::locale(),
  //         boost::property_tree::xml_writer_make_settings<std::string>(' ',
  //         4));
};

BOOST_AUTO_TEST_CASE(AmpTreeCorrespondence) {
    boost::property_tree::ptree tr;
  std::stringstream modelStream;
  // Construct particle list from XML tree
  modelStream << HelicityTestParticles;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, tr);

  modelStream.clear();
  tr = boost::property_tree::ptree();
  // Construct Kinematics from XML tree
  modelStream << HelicityTestKinematics;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto kin = std::make_shared<HelicityKinematics>(
      partL, tr.get_child("HelicityKinematics"));

  modelStream.clear();
  tr = boost::property_tree::ptree();
  modelStream << HelicityTestModel;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  // Create amplitude from property_tree
  auto intens = ComPWA::Physics::IncoherentIntensity::Factory(
      partL, kin, tr.get_child("IncoherentIntensity"));

  ParameterList list;
  intens->parameters(list);
  LOG(info) << "List of parameters: ";
  for( auto p : list.doubleParameters() )
    LOG(info) << p->to_str();
  BOOST_CHECK_EQUAL(list.doubleParameters().size(), 14);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      partL, kin->initialState(), kin->finalState(), 123));
  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::Data());
  
  ComPWA::Tools::GeneratePhsp(200, gen, sample);

  auto phspSample =
      std::make_shared<std::vector<DataPoint>>(sample->GetDataPoints(kin));
  intens->setPhspSample(phspSample, phspSample);

  LOG(info) << "Loop over phsp events....";
  for (auto i : sample->GetEvents()) {
    ComPWA::DataPoint p;
    try {
      kin->convert(i, p);
    } catch (std::exception &ex) {
      // Test if events outside the phase space boundaries are generated
      LOG(error) << "Event outside phase space. This should not happen since "
                    "we use a Monte-Carlo sample!";
      BOOST_FAIL("Event outside phase space. This should not happen since "
                 "we use a Monte-Carlo sample!");
      continue;
    }

    double w = intens->intensity(p);
    i.setWeight(w);
    LOG(info) << "point = " << p << " intensity = " << w;
  }

  ComPWA::ParameterList sampleList(sample->GetListOfData(kin));
  // Testing function tree
  auto tree = intens->tree(kin, sampleList, sampleList, sampleList,
                           kin->numVariables());
  tree->parameter();

  // TODO: implement checks to ensure that amplitude calculation by FunctionTree
  //      and by Evaluate() are the same

  std::stringstream printTree;
  printTree << tree;
  LOG(info) << std::endl << printTree.str();
};
BOOST_AUTO_TEST_CASE(AmpMinimizationUsingTree) {}
BOOST_AUTO_TEST_CASE(AmpMinimization) {}

BOOST_AUTO_TEST_SUITE_END()
