
        
  
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
#include "Tools/RunManager.hpp"
#include "Tools/RootGenerator.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

using namespace ComPWA::Physics::HelicityFormalism;

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <iostream>

using namespace ComPWA;

BOOST_AUTO_TEST_SUITE(HelicityFormalism)

BOOST_AUTO_TEST_CASE(KinematicsConstructionFromXML) {
  ComPWA::Logging log("", boost::log::trivial::severity_level::trace);

  // Construct HelicityKinematics from XML tree
  boost::property_tree::ptree tr;
  boost::property_tree::xml_parser::read_xml(
      "../HelicityFormalismTest-input.xml", tr);

  ComPWA::PhysConst::CreateInstance(tr);
  HelicityKinematics::CreateInstance(tr.get_child("HelicityKinematics"));

  auto inst = HelicityKinematics::Instance();
  BOOST_CHECK_EQUAL(inst->GetPhspVolume(), 0.123);
  BOOST_CHECK_EQUAL(inst->GetFinalState().size(), 3);
  BOOST_CHECK_EQUAL(inst->GetInitialState().size(), 1);
}

BOOST_AUTO_TEST_CASE(ConstructionFromXML) {
  boost::property_tree::ptree tr;
  boost::property_tree::xml_parser::read_xml(
      "../HelicityFormalismTest-input.xml", tr);

  // Due to the structure of Boost.UnitTest the instances already exist from
  // previous test
  //  ComPWA::Logging log("", boost::log::trivial::severity_level::trace);
  //  ComPWA::PhysConst::CreateInstance(tr);
  //  HelicityKinematics::CreateInstance(tr.get_child("HelicityKinematics"));

  // Create amplitude from property_tree
  auto intens =
      IncoherentIntensity::Factory(tr.get_child("IncoherentIntensity"));

  // Save amplitude to property_tree
  boost::property_tree::ptree ptout;
  ptout.add_child("IncoherentIntensity", IncoherentIntensity::Save(intens));

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
  boost::property_tree::xml_parser::write_xml(
      "../HelicityFormalismTest-output.xml", ptout, std::locale(),
      boost::property_tree::xml_writer_make_settings<std::string>(' ', 4,
                                                                  "utf-8"));
};

BOOST_AUTO_TEST_CASE(AmpTreeCorrespondence) {
  boost::property_tree::ptree tr;
  boost::property_tree::xml_parser::read_xml(
      "../HelicityFormalismTest-input.xml", tr);

  // Due to the structure of Boost.UnitTest the instances already exist from
  // previous test
  //  ComPWA::Logging log("", boost::log::trivial::severity_level::trace);
  //  ComPWA::PhysConst::CreateInstance(tr);
  //  HelicityKinematics::CreateInstance(tr.get_child("HelicityKinematics"));

  // Create amplitude
  auto intens =
      IncoherentIntensity::Factory(tr.get_child("IncoherentIntensity"));

  ParameterList list;
  intens->GetParameters(list);
  LOG(info) << "List of parameters: " << std::endl << list;
  BOOST_CHECK_EQUAL(list.GetNDouble(), 14);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(123));
  std::shared_ptr<ComPWA::DataReader::Data> sample(
      new ComPWA::DataReader::RootReader());

  ComPWA::RunManager r;
  r.SetGenerator(gen);
  r.SetPhspSample(sample);
  r.GeneratePhsp(200);

  auto phspSample = std::make_shared<std::vector<dataPoint>>(
      r.GetPhspSample()->GetDataPoints());
  intens->SetPhspSample(phspSample, phspSample);

  LOG(info) << "Loop over phsp events....";
  for (auto i : sample->GetEvents()) {
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
    i.SetWeight(w);
    LOG(info) << "point = " << p << " intensity = " << w;
  }

  ComPWA::ParameterList sampleList(sample->GetListOfData());
  // Testing function tree
  auto tree = intens->GetTree(sampleList, sampleList, sampleList);
  tree->recalculate();

  // TODO: implement checks to ensure that amplitude calculation by FunctionTree
  //      and by Evaluate() are the same

  std::stringstream printTree;
  printTree << tree;
  LOG(info) << std::endl << printTree.str();
};
BOOST_AUTO_TEST_CASE(AmpMinimizationUsingTree) {}
BOOST_AUTO_TEST_CASE(AmpMinimization) {}

BOOST_AUTO_TEST_SUITE_END()
