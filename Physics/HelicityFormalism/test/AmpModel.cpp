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
#include "Physics/DecayDynamics/RelativisticBreitWigner.hpp"

// Amplitude models
#include "Physics/HelicityFormalism/test/AmpModelTest.hpp"

using namespace ComPWA;
using namespace ComPWA::Physics::HelicityFormalism;

BOOST_AUTO_TEST_SUITE(HelicityFormalism)

BOOST_AUTO_TEST_CASE(KinematicsConstructionFromXML) {
  ComPWA::Logging log("", "trace");

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
  BOOST_CHECK_EQUAL(kin->getKinematicsProperties().FinalState.size(), 3);
  BOOST_CHECK_EQUAL(kin->getKinematicsProperties().InitialState.size(), 1);
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
  auto intens = std::make_shared<ComPWA::Physics::IncoherentIntensity>(
      partL, kin, tr.get_child("Intensity"));

  // Save amplitude to property_tree
  boost::property_tree::ptree ptout;
  ptout.add_child("Intensity",intens->save());

  //  if (ptout != tr) {
  //    BOOST_CHECK(false);
  //    LOG(ERROR)
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

BOOST_AUTO_TEST_CASE(PartialAmplitudeTreeConcordance) {
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
  modelStream << PartialAmplitudeTestModel;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto helDecay = std::make_shared<
      ComPWA::Physics::HelicityFormalism::HelicityDecay>(
      partL, kin, tr.get_child("PartialAmplitude"));

  ParameterList list;
  helDecay->parameters(list);
  LOG(INFO) << "List of parameters: ";
  for( auto p : list.doubleParameters() )
    LOG(INFO) << p->to_str();
  BOOST_CHECK_EQUAL(list.doubleParameters().size(), 5);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      partL, kin->getKinematicsProperties().InitialState,
      kin->getKinematicsProperties().FinalState, 123));
  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::Data());

  std::shared_ptr<ComPWA::DataReader::Data> toySample(
      new ComPWA::DataReader::Data());
  
  ComPWA::Tools::generatePhsp(20, gen, sample);
  ComPWA::Tools::generatePhsp(20000, gen, toySample);

  auto phspSample =
      std::make_shared<std::vector<DataPoint>>(sample->dataPoints(kin));
  auto toyPhspSample =
      std::make_shared<std::vector<DataPoint>>(toySample->dataPoints(kin));
  helDecay->setPhspSample(toyPhspSample);

  ComPWA::ParameterList sampleList(sample->dataList(kin));
  ComPWA::ParameterList toySampleList(toySample->dataList(kin));
  // Testing function tree
  auto tree = helDecay->tree(kin, sampleList, toySampleList);
  tree->parameter();
  LOG(INFO)<<tree->print();
  
  LOG(INFO) << "Loop over phsp events....";
  for (size_t i = 0; i < sample->numEvents(); i++) {
    ComPWA::DataPoint p;
    try {
      kin->convert(sample->event(i), p);
    } catch (std::exception &ex) {
      // Test if events outside the phase space boundaries are generated
      LOG(ERROR) << "Event outside phase space. This should not happen since "
                    "we use a Monte-Carlo sample!";
      BOOST_FAIL("Event outside phase space. This should not happen since "
                 "we use a Monte-Carlo sample!");
      continue;
    }
    // Intensity without function tree
    auto intensityNoTree = helDecay->evaluate(p);
    
    // Intensity calculated using function tree
    auto tmp = tree->head()->parameter();
    auto intensitiesTree =
        std::dynamic_pointer_cast<Value<std::vector<std::complex<double>>>>(tmp);
    std::complex<double> intensityTree = intensitiesTree->values().at(i);

    //// Not available in boost 1.54
    //BOOST_TEST(intensityNoTree.real() == intensityTree.real(),
    //            boost::test_tools::tolerance(0.000001));
    //BOOST_TEST(intensityNoTree.imag() == intensityTree.imag(),
    //            boost::test_tools::tolerance(0.000001));
    BOOST_CHECK_CLOSE(intensityNoTree.real(), intensityTree.real(), 0.0000001);
    BOOST_CHECK_CLOSE(intensityNoTree.imag(), intensityTree.imag(), 0.0000001);

    LOG(DEBUG) << "point = " << p << " intensity = " << intensityNoTree
              << " intensity tree = " << intensityTree;
  }
};

BOOST_AUTO_TEST_CASE(RelBWTreeConcordance) {
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

  // we test the breit-wigner but we need to construct the full amplitude in
  // order to register the variables in HelicityKinematics
  modelStream.clear();
  tr = boost::property_tree::ptree();
  modelStream << HelicityTestModel;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  // Create amplitude from property_tree
  auto intens = std::make_shared<ComPWA::Physics::IncoherentIntensity>(
      partL, kin, tr.get_child("Intensity"));

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      partL, kin->getKinematicsProperties().InitialState,
      kin->getKinematicsProperties().FinalState, 123));
  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::Data());

  ComPWA::Tools::generatePhsp(20, gen, sample);

  auto relBW =
      std::make_shared<ComPWA::Physics::DecayDynamics::RelativisticBreitWigner>(
          "omega", std::make_pair("pi0", "gamma"), partL);
  ComPWA::ParameterList sampleList(sample->dataList(kin));
  // Testing function tree
  auto tree = relBW->tree(sampleList, 3);
  tree->parameter();
  LOG(INFO)<<tree->print();
  
  LOG(INFO) << "Loop over phsp events....";
  for (size_t i = 0; i < sample->numEvents(); i++) {
    ComPWA::DataPoint p;
    try {
      kin->convert(sample->event(i), p);
    } catch (std::exception &ex) {
      // Test if events outside the phase space boundaries are generated
      LOG(ERROR) << "Event outside phase space. This should not happen since "
                    "we use a Monte-Carlo sample!";
      BOOST_FAIL("Event outside phase space. This should not happen since "
                 "we use a Monte-Carlo sample!");
      continue;
    }
    // Intensity without function tree
    auto intensityNoTree = relBW->evaluate(p, 3);
    
    // Intensity calculated using function tree
    auto tmp = tree->head()->parameter();
    auto intensitiesTree =
        std::dynamic_pointer_cast<Value<std::vector<std::complex<double>>>>(tmp);
    std::complex<double> intensityTree = intensitiesTree->values().at(i);
    
    BOOST_CHECK_EQUAL(intensityNoTree, intensityTree);
    LOG(DEBUG) << "point = " << p << " intensity = " << intensityNoTree
              << " intensity tree = " << intensityTree;
  }
};

BOOST_AUTO_TEST_CASE(IncoherentTreeConcordance) {
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
  auto intens = std::make_shared<ComPWA::Physics::IncoherentIntensity>(
      partL, kin, tr.get_child("Intensity"));

  ParameterList list;
  intens->parameters(list);
  LOG(INFO) << "List of parameters: ";
  for( auto p : list.doubleParameters() )
    LOG(INFO) << p->to_str();
  BOOST_CHECK_EQUAL(list.doubleParameters().size(), 14);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      partL, kin->getKinematicsProperties().InitialState,
      kin->getKinematicsProperties().FinalState, 123));
  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::Data());

  std::shared_ptr<ComPWA::DataReader::Data> toySample(
      new ComPWA::DataReader::Data());
  
  ComPWA::Tools::generatePhsp(20, gen, sample);
  ComPWA::Tools::generatePhsp(20000, gen, toySample);

  auto vToySample =
      std::make_shared<std::vector<DataPoint>>(toySample->dataPoints(kin));
  intens->setPhspSample(vToySample, vToySample);

  ComPWA::ParameterList sampleList(sample->dataList(kin));
  ComPWA::ParameterList toySampleList(toySample->dataList(kin));
  // Testing function tree
  auto tree = intens->tree(kin, sampleList, toySampleList, toySampleList,
                           kin->numVariables());
  tree->parameter();
  LOG(INFO)<<tree->print();
  
  LOG(INFO) << "Loop over phsp events....";
  for (size_t i = 0; i < sample->numEvents(); i++) {
    ComPWA::DataPoint p;
    try {
      kin->convert(sample->event(i), p);
    } catch (std::exception &ex) {
      // Test if events outside the phase space boundaries are generated
      LOG(ERROR) << "Event outside phase space. This should not happen since "
                    "we use a Monte-Carlo sample!";
      BOOST_FAIL("Event outside phase space. This should not happen since "
                 "we use a Monte-Carlo sample!");
      continue;
    }
    // Intensity without function tree
    auto intensityNoTree = intens->intensity(p);
    
    // Intensity calculated using function tree
    auto tmp = tree->head()->parameter();
    auto intensitiesTree =
        std::dynamic_pointer_cast<Value<std::vector<double>>>(tmp);
    double intensityTree = intensitiesTree->values().at(i);
    
    //// Not available in boost 1.54
    //BOOST_TEST(intensityNoTree == intensityTree,
    //           boost::test_tools::tolerance(0.000001));
    BOOST_CHECK_CLOSE(intensityNoTree, intensityTree, 0.0000001);
    LOG(INFO) << "point = " << p << " intensity = " << intensityNoTree
              << " intensity tree = " << intensityTree;
  }
};

BOOST_AUTO_TEST_CASE(SeqPartialAmplitudeTreeConcordance) {
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
  modelStream << SeqPartialAmplitudeTestModel;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto seqAmp = std::make_shared<
      ComPWA::Physics::HelicityFormalism::SequentialPartialAmplitude>(
      partL, kin, tr.get_child("Amplitude"));

  ParameterList list;
  seqAmp->parameters(list);
  LOG(INFO) << "List of parameters: ";
  for( auto p : list.doubleParameters() )
    LOG(INFO) << p->to_str();
  BOOST_CHECK_EQUAL(list.doubleParameters().size(), 12);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      partL, kin->getKinematicsProperties().InitialState,
      kin->getKinematicsProperties().FinalState, 123));
  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::Data());

  std::shared_ptr<ComPWA::DataReader::Data> toySample(
      new ComPWA::DataReader::Data());
  
  ComPWA::Tools::generatePhsp(20, gen, sample);
  ComPWA::Tools::generatePhsp(20000, gen, toySample);

  auto phspSample =
      std::make_shared<std::vector<DataPoint>>(sample->dataPoints(kin));
  auto toyPhspSample =
      std::make_shared<std::vector<DataPoint>>(toySample->dataPoints(kin));
  seqAmp->setPhspSample(toyPhspSample);

  ComPWA::ParameterList sampleList(sample->dataList(kin));
  ComPWA::ParameterList toySampleList(toySample->dataList(kin));
  // Testing function tree
  auto tree = seqAmp->tree(kin, sampleList, toySampleList);
  tree->parameter();
  LOG(INFO)<<tree->print();
  
  LOG(INFO) << "Loop over phsp events....";
  for (size_t i = 0; i < sample->numEvents(); i++) {
    ComPWA::DataPoint p;
    try {
      kin->convert(sample->event(i), p);
    } catch (std::exception &ex) {
      // Test if events outside the phase space boundaries are generated
      LOG(ERROR) << "Event outside phase space. This should not happen since "
                    "we use a Monte-Carlo sample!";
      BOOST_FAIL("Event outside phase space. This should not happen since "
                 "we use a Monte-Carlo sample!");
      continue;
    }
    // Intensity without function tree
    auto intensityNoTree = seqAmp->evaluate(p);
    
    // Intensity calculated using function tree
    auto tmp = tree->head()->parameter();
    auto intensitiesTree =
        std::dynamic_pointer_cast<Value<std::vector<std::complex<double>>>>(tmp);
    std::complex<double> intensityTree = intensitiesTree->values().at(i);
    
    //// Not available in boost 1.54
    //BOOST_TEST(intensityNoTree.real() == intensityTree.real(),
    //           boost::test_tools::tolerance(0.000001));
    //BOOST_TEST(intensityNoTree.imag() == intensityTree.imag(),
    //           boost::test_tools::tolerance(0.000001));
    BOOST_CHECK_CLOSE(intensityNoTree.real(), intensityTree.real(), 0.0000001);
    BOOST_CHECK_CLOSE(intensityNoTree.imag(), intensityTree.imag(), 0.0000001);

    LOG(INFO) << "point = " << p << " intensity = " << intensityNoTree
              << " intensity tree = " << intensityTree;
  }
};



BOOST_AUTO_TEST_SUITE_END()
