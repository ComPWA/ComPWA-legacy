// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE UpdatePTreeParameterTest

#include <cmath>
#include <vector>

#include "Core/Logging.hpp"

#include <boost/foreach.hpp>
#include <boost/locale/utf.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

#include "Tools/UpdatePTreeParameter.hpp"

BOOST_AUTO_TEST_SUITE(ToolsTest)

const std::string JpsiDecayTree = R"####(
<Intensity Class='StrengthIntensity' Name='incoherent_with_strength'>
  <Parameter Class='Double' Type='Strength' Name='strength_incoherent'>
    <Value>1</Value>
    <Fix>true</Fix>
  </Parameter>
  <Intensity Class='IncoherentIntensity' Name='incoherent'>
    <Intensity Class='CoherentIntensity' Name='coherent_0'>
      <Amplitude Class='CoefficientAmplitude' Name='jpsi_1_to_gamma_1+f0_0_L_0.0_S_1.0'>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_0.0_S_1.0'>
          <Value>1.0</Value>
          <Fix>false</Fix>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_0.0_S_1.0'>
          <Value>0.0</Value>
          <Fix>false</Fix>
        </Parameter>
      </Amplitude>
      <Amplitude Class='CoefficientAmplitude' Name='jpsi_1_to_gamma_1+f0_0_L_2.0_S_1.0'>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_2.0_S_1.0'>
          <Value>1.0</Value>
          <Fix>false</Fix>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_2.0_S_1.0'>
          <Value>0.0</Value>
          <Fix>false</Fix>
        </Parameter>
      </Amplitude>
    </Intensity>

    <Intensity Class='CoherentIntensity' Name='coherent_1'>
      <Amplitude Class='CoefficientAmplitude' Name='jpsi_-1_to_gamma_1+f0_0_L_0.0_S_1.0'>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_0.0_S_1.0'>
          <Value>1.0</Value>
          <Fix>false</Fix>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_0.0_S_1.0'>
          <Value>0.0</Value>
          <Fix>false</Fix>
        </Parameter>
      </Amplitude>
      <Amplitude Class='CoefficientAmplitude' Name='jpsi_-1_to_gamma_1+f0_0_L_2.0_S_1.0'>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_2.0_S_1.0'>
          <Value>1.0</Value>
          <Fix>false</Fix>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_2.0_S_1.0'>
          <Value>0.0</Value>
          <Fix>false</Fix>
        </Parameter>
      </Amplitude>
    </Intensity>

  </Intensity>
</Intensity>
)####";

BOOST_AUTO_TEST_CASE(UpdatePTreeParameterTest) {
  ComPWA::Logging Log("trace", "");

  LOG(INFO) << "Now check functions in "
            << "ComPWA::Tools/UpdatePTreeParameter.cpp ...";

  boost::property_tree::ptree ModelTree;
  std::stringstream ModelStream;

  // Model intensity
  ModelStream.clear();
  ModelTree = boost::property_tree::ptree();
  ModelStream << JpsiDecayTree;
  boost::property_tree::xml_parser::read_xml(ModelStream, ModelTree);

  // there are five unique parameters
  // type Strength, number = 1
  std::string Strength("strength_incoherent");
  // type Magnitude, number = 2
  std::string MagnitudeL2S1("Magnitude_jpsi_to_gamma+f0_L_2.0_S_1.0");
  // type Phase, number = 2
  std::string PhaseL2S1("Phase_jpsi_to_gamma+f0_L_2.0_S_1.0");
  // type Magnitude, number = 2
  std::string MagnitudeL0S1("Magnitude_jpsi_to_gamma+f0_L_0.0_S_1.0");
  // type Phase, number = 2
  std::string PhaseL0S1("Phase_jpsi_to_gamma+f0_L_0.0_S_1.0");

  boost::property_tree::ptree RootTree(ModelTree.get_child("Intensity"));

  void check(boost::property_tree::ptree & Tree, const std::string &KeyType,
             const std::string &KeyValue, int &Count, bool CheckValue,
             double Value, bool CheckBounds, double Min, double Max,
             bool CheckFix, bool Fix);

  int Count = 0;

  ComPWA::Tools::updateParameterRangeByType(RootTree, "Magnitude", 0, 10);

  check(RootTree, "Type", "Magnitude", Count,
        //   dummy args  real  args   dummy args
        false, 999, true, 0, 10, false, false);
  LOG(INFO) << "4 Magnitude Parameters range are changed to [0, 10], " << Count
            << " Magnitude Parameters' range are right changed.";
  BOOST_CHECK_EQUAL(Count, 4);

  ComPWA::Tools::updateParameterRangeByName(RootTree, Strength, 0, 1);
  double MinTemp = RootTree.get<double>("Parameter.Min");
  double MaxTemp = RootTree.get<double>("Parameter.Max");
  LOG(INFO) << "Range of " << Strength << " is set to [0, 1], ";
  LOG(INFO) << "Now the range is [" << MinTemp << ", " << MaxTemp << "].";
  BOOST_CHECK_EQUAL(MinTemp, 0);
  BOOST_CHECK_EQUAL(MaxTemp, 1);

  ComPWA::Tools::updateParameterValue(RootTree, Strength, 10);
  double ValTemp = RootTree.get<double>("Parameter.Value");
  LOG(INFO) << "valTemp = " << ValTemp;
  BOOST_CHECK_EQUAL(ValTemp, 10);
  LOG(INFO) << "Value of " << Strength << " is set to 10, current value is "
            << ValTemp;

  ComPWA::Tools::fixParameter(RootTree, PhaseL0S1);
  Count = 0;
  check(RootTree, "Name", PhaseL0S1, Count, false, 999, false, 999, 999, true,
        true);
  LOG(INFO) << "2 parameters " << PhaseL0S1 << " are fixed, current number of "
            << "fixed " << PhaseL0S1 << " is " << Count;
  BOOST_CHECK_EQUAL(2, Count);

  ComPWA::Tools::releaseParameter(RootTree, PhaseL0S1);
  Count = 0;
  check(RootTree, "Name", PhaseL0S1, Count, false, 999, false, 999, 999, true,
        true);
  LOG(INFO) << "Now after realeasParameter, there should " << Count << " fixed "
            << PhaseL0S1 << " parameter";
  BOOST_CHECK_EQUAL(0, Count);

  auto FitStrength =
      ComPWA::FitParameter<double>(Strength, 0.5, 0.0, 1.0, false);
  auto FitMagnitudeL0S1 =
      ComPWA::FitParameter<double>(MagnitudeL0S1, 1.0, 0.0, 10.0, false);
  auto FitPhaseL0S1 =
      ComPWA::FitParameter<double>(PhaseL0S1, 0.0, -3.14, 3.14, false);
  auto FitMagnitudeL2S1 =
      ComPWA::FitParameter<double>(MagnitudeL2S1, 3.0, 0.0, 10.0, false);
  auto FitPhaseL2S1 =
      ComPWA::FitParameter<double>(PhaseL2S1, 1.2, -3.14, 3.14, false);
  std::vector<ComPWA::FitParameter<double>> FitParameters(
      {FitStrength, FitMagnitudeL0S1, FitPhaseL0S1, FitMagnitudeL2S1,
       FitPhaseL2S1});

  // RootTree = ModelTree.get_child("Intensity");
  ComPWA::Tools::updateParameter(RootTree, FitParameters);
  Count = 0;
  check(RootTree, "Name", Strength, Count, true, 0.5, true, 0.0, 1.0, true,
        false);
  BOOST_CHECK_EQUAL(Count, 1);
  Count = 0;
  check(RootTree, "Name", MagnitudeL0S1, Count, true, 1.0, true, 0.0, 10.0,
        true, false);
  BOOST_CHECK_EQUAL(Count, 2);
  Count = 0;
  check(RootTree, "Name", PhaseL0S1, Count, true, 0.0, true, -3.14, 3.14, true,
        false);
  BOOST_CHECK_EQUAL(Count, 2);

  Count = 0;
  check(RootTree, "Name", MagnitudeL2S1, Count, true, 3.0, true, 0.0, 10.0,
        true, false);
  BOOST_CHECK_EQUAL(Count, 2);
  Count = 0;
  check(RootTree, "Name", PhaseL2S1, Count, true, 1.2, true, -3.14, 3.14, true,
        false);
  BOOST_CHECK_EQUAL(Count, 2);
}
void check(boost::property_tree::ptree &Tree, const std::string &KeyType,
           const std::string &KeyValue, int &Count, bool CheckValue,
           double Value, bool CheckBounds, double Min, double Max,
           bool CheckFix, bool Fix) {

  for (const auto &Node : Tree) {
    // if this node is not parameter, recursively check this node
    if (Node.first != "Parameter") {
      auto Child = Node.second;
      check(Child, KeyType, KeyValue, Count, CheckValue, Value, CheckBounds,
            Min, Max, CheckFix, Fix);
      continue;
    }

    // if this is a parameter, but not the target one, continue
    if (KeyValue != Node.second.get<std::string>("<xmlattr>." + KeyType)) {
      continue;
    }

    bool ValueOK(true), BoundsOK(true), FixOK(true);
    // check the target parameter's properties and Count the right one
    if (CheckValue && abs(Value - Node.second.get<double>("Value")) > 1e-6) {
      ValueOK = false;
    }
    if (CheckBounds && (abs(Min - Node.second.get<double>("Min")) > 1e-6 ||
                        abs(Max - Node.second.get<double>("Max")) > 1e-6)) {
      BoundsOK = false;
    }
    if (CheckFix && Node.second.get<bool>("Fix") != Fix) {
      FixOK = false;
    }

    if (ValueOK && BoundsOK && FixOK)
      Count++;
  }
}

BOOST_AUTO_TEST_SUITE_END()
