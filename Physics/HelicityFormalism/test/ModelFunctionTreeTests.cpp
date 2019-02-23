// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// This can only be define once within the same library ?!
#define BOOST_TEST_MODULE HelicityFormalism

#include <vector>

#include "Core/Intensity.hpp"
#include "Core/Logging.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Physics/Amplitude.hpp"
#include "Physics/Dynamics/RelativisticBreitWigner.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Tools/Generate.hpp"
#include "Tools/RootGenerator.hpp"

#include <boost/foreach.hpp>
#include <boost/locale/utf.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

using namespace ComPWA;
using namespace ComPWA::Physics::HelicityFormalism;

BOOST_AUTO_TEST_SUITE(HelicityFormalism)

const std::string HelicityTestParticles = R"####(
<ParticleList>
  <Particle Name='pi0'>
    <Pid>111</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_pi0'>
      <Value>0.1349766</Value>
      <Error>0.0000006</Error>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='1'/>
  </Particle>
  <Particle Name='gamma'>
    <Pid>22</Pid>
    <Parameter Class='Double' Type='Mass' Name='mass_gamma'>
      <Value>0.</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
  </Particle>
  <Particle Name='jpsi'>
    <Pid>443</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_jpsi'>
      <Value>3.0969</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='0' />
      <Parameter Class='Double' Type='Width' Name='Width_jpsi'>
        <Value>0.0000929</Value>
        <Error>0.0000028</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_jpsi'>
        <Value>2.5</Value>
        <Fix>true</Fix>
        <Min>2.0</Min>
        <Max>3.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name='omega'>
    <Pid>223</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_omega'>
      <Value>0.78265</Value>
      <Fix>true</Fix>
      <Min>0.5</Min>
      <Max>1.5</Max>
      <Error>0.00012</Error>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='0' />
      <Parameter Class='Double' Type='Width' Name='Width_omega'>
        <Value>0.01849</Value>
        <Fix>true</Fix>
        <Min>0.0</Min>
        <Max>1.0</Max>
        <Error>0.00008</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_omega'>
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
        <Error>0</Error>
      </Parameter>
    </DecayInfo>
  </Particle>
</ParticleList>
)####";

const std::string HelicityTestKinematics = R"####(
<HelicityKinematics>
  <PhspVolume>0.123</PhspVolume>
  <InitialState>
    <Particle Name='jpsi' PositionIndex='0'/>
  </InitialState>
  <FinalState>
    <Particle Name='pi0' Id='0'/>
    <Particle Name='gamma' Id='1'/>
    <Particle Name='pi0' Id='2'/>
  </FinalState>
</HelicityKinematics>
)####";

const std::string HelicityTestModel = R"####(
<Intensity Class='NormalizedIntensity' Name='jpsiToPi0Pi0Gamma_norm'>
  <Intensity Class='StrengthIntensity' Name='jpsiToPi0Pi0Gamma_strength'>
  <Parameter Class='Double' Type='Strength' Name='Strength_jpsiToPi0Pi0Gamma'>
    <Value>0.99</Value>
    <Fix>true</Fix>
  </Parameter>
  <Intensity Class='CoherentIntensity' Name='jpsiToPi0Pi0Gamma'>
    <Amplitude Class='CoefficientAmplitude' Name='omega_coeff'>
      <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsitoOmegaPi0'>
        <Value>1.0</Value>
        <Fix>true</Fix>
        <Min>0.5</Min>
        <Max>1.5</Max>
        <Error>0</Error>
      </Parameter>
      <Parameter Class='Double' Type='Phase' Name='Phase_jpsitoOmegaPi0'>
        <Value>1.0</Value>
        <Fix>true</Fix>
        <Min>0.5</Min>
        <Max>1.5</Max>
        <Error>0</Error>
      </Parameter>
      <Amplitude Class='SequentialAmplitude' Name='omega'>
        <Amplitude Class="HelicityDecay" Name='jpsitoOmegaPi0'>
          <DecayParticle Name='jpsi' Helicity='+1' />
          <DecayProducts>
            <Particle Name='omega' FinalState='0 1' Helicity='+1' />
            <Particle Name='pi0' FinalState='2' Helicity='0' />
          </DecayProducts>
        </Amplitude>
        <Amplitude Class="HelicityDecay" Name="omegatoPi0G">
          <DecayParticle Name='omega' Helicity='+1' />
          <RecoilSystem FinalState='2' />
          <DecayProducts>
            <Particle Name='gamma' FinalState='1' Helicity='+1' />
            <Particle Name='pi0' FinalState='0' Helicity='0' />
          </DecayProducts>
        </Amplitude>
      </Amplitude>
    </Amplitude>
  </Intensity>
  </Intensity>
</Intensity>
)####";

const std::string SequentialAmplitudeTestModel = R"####(
    <Amplitude Class='CoefficientAmplitude' Name='omega_norm'>
      <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsitoOmegaPi0'>
        <Value>1.0</Value>
        <Fix>true</Fix>
        <Min>0.5</Min>
        <Max>1.5</Max>
        <Error>0</Error>
      </Parameter>
      <Parameter Class='Double' Type='Phase' Name='Phase_jpsitoOmegaPi0'>
        <Value>1.0</Value>
        <Fix>true</Fix>
        <Min>0.5</Min>
        <Max>1.5</Max>
        <Error>0</Error>
      </Parameter>
      <Amplitude Class='SequentialAmplitude' Name='omega'>
        <Amplitude Class="HelicityDecay" Name='jpsitoOmegaPi0'>
          <DecayParticle Name='jpsi' Helicity='+1' />
          <DecayProducts>
            <Particle Name='omega' FinalState='0 1' Helicity='+1' />
            <Particle Name='pi0' FinalState='2' Helicity='0' />
          </DecayProducts>
        </Amplitude>
        <Amplitude Class="HelicityDecay" Name="omegatoPi0G">
          <DecayParticle Name='omega' Helicity='+1' />
          <RecoilSystem FinalState='2' />
          <DecayProducts>
            <Particle Name='gamma' FinalState='1' Helicity='+1' />
            <Particle Name='pi0' FinalState='0' Helicity='0' />
          </DecayProducts>
        </Amplitude>
      </Amplitude>
    </Amplitude>
)####";

const std::string HelicityDecayTestModel = R"####(
    <Amplitude Class='CoefficientAmplitude' Name='omega_norm'>
      <Parameter Class='Double' Type='Magnitude' Name='Magnitude_omegaToPi0Gamma'>
        <Value>1.0</Value>
        <Fix>true</Fix>
        <Min>0.5</Min>
        <Max>1.5</Max>
        <Error>0</Error>
      </Parameter>
      <Parameter Class='Double' Type='Phase' Name='Phase_omegaToPi0Gamma'>
        <Value>1.0</Value>
        <Fix>true</Fix>
        <Min>0.5</Min>
        <Max>1.5</Max>
        <Error>0</Error>
      </Parameter>
      <Amplitude Class="HelicityDecay" Name="omegatoPi0G">
        <DecayParticle Name='omega' Helicity='+1' />
        <RecoilSystem FinalState='2' />
        <DecayProducts>
          <Particle Name='gamma' FinalState='1' Helicity='+1' />
          <Particle Name='pi0' FinalState='0' Helicity='0' />
        </DecayProducts>
      </Amplitude>
    </Amplitude>
)####";

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

  ComPWA::Physics::IntensityBuilderXML Builder;
  auto kin = Builder.createHelicityKinematics(
      partL, tr.get_child("HelicityKinematics"));

  BOOST_CHECK_EQUAL(kin->phspVolume(), 0.123);
  BOOST_CHECK_EQUAL(kin->getParticleStateTransitionKinematicsInfo()
                        .getFinalStateParticleCount(),
                    3);
  BOOST_CHECK_EQUAL(kin->getParticleStateTransitionKinematicsInfo()
                        .getInitialStateInvariantMass(),
                    partL->at("jpsi").GetMass());
}

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

  ComPWA::Physics::IntensityBuilderXML Builder;
  auto kin = Builder.createHelicityKinematics(
      partL, tr.get_child("HelicityKinematics"));

  modelStream.clear();
  tr = boost::property_tree::ptree();
  modelStream << HelicityDecayTestModel;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);

  std::shared_ptr<ComPWA::Physics::NamedAmplitude> helDecay =
      Builder.createAmplitude(partL, kin, tr.get_child("Amplitude"));

  ParameterList list;
  helDecay->addUniqueParametersTo(list);
  LOG(INFO) << "List of parameters: ";
  for (auto p : list.doubleParameters())
    LOG(INFO) << p->to_str();
  BOOST_CHECK_EQUAL(list.doubleParameters().size(), 5);

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      kin->getParticleStateTransitionKinematicsInfo(), 123));
  std::shared_ptr<ComPWA::Data::DataSet> sample(
      ComPWA::Tools::generatePhsp(20, gen));

  sample->convertEventsToParameterList(kin);

  // Testing function tree
  auto tree = helDecay->createFunctionTree(sample->getParameterList(), "");
  auto tmp = tree->parameter();
  LOG(INFO) << tree->print();
  // Intensity calculated using function tree
  auto intensitiesTree =
      std::dynamic_pointer_cast<Value<std::vector<std::complex<double>>>>(tmp);

  unsigned int counter(0);
  LOG(INFO) << "Loop over events....";
  for (auto const &x : sample->getDataPointList()) {
    // Intensity without function tree
    auto intensityNoTree = helDecay->evaluate(x);

    std::complex<double> intensityTree =
        intensitiesTree->values().at(counter++);

    //// Not available in boost 1.54
    // BOOST_TEST(intensityNoTree.real() == intensityTree.real(),
    //            boost::test_tools::tolerance(0.000001));
    // BOOST_TEST(intensityNoTree.imag() == intensityTree.imag(),
    //            boost::test_tools::tolerance(0.000001));
    BOOST_CHECK_CLOSE(intensityNoTree.real(), intensityTree.real(), 0.0000001);
    BOOST_CHECK_CLOSE(intensityNoTree.imag(), intensityTree.imag(), 0.0000001);

    LOG(DEBUG) << "point = " << x << " intensity = " << intensityNoTree
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
  ComPWA::Physics::IntensityBuilderXML Builder;
  auto kin = Builder.createHelicityKinematics(
      partL, tr.get_child("HelicityKinematics"));

  modelStream.clear();
  tr = boost::property_tree::ptree();
  modelStream << HelicityDecayTestModel;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);

  // Generate sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      kin->getParticleStateTransitionKinematicsInfo(), 123));
  std::shared_ptr<ComPWA::Data::DataSet> sample(
      ComPWA::Tools::generatePhsp(20, gen));

  kin->addSubSystem({0}, {1}, {2}, {});
  sample->convertEventsToParameterList(kin);

  auto relBW =
      std::make_shared<ComPWA::Physics::Dynamics::RelativisticBreitWigner>(
          "omega", std::make_pair("pi0", "gamma"), partL);

  // Testing function tree
  auto tree = relBW->createFunctionTree(sample->getParameterList(), 0, "");
  auto tmp = tree->parameter();
  LOG(INFO) << tree->print();
  // Intensity calculated using function tree
  auto intensitiesTree =
      std::dynamic_pointer_cast<Value<std::vector<std::complex<double>>>>(tmp);

  unsigned int counter(0);
  LOG(INFO) << "Loop over events....";
  for (auto const &x : sample->getDataPointList()) {
    // Intensity without function tree
    auto intensityNoTree = relBW->evaluate(x, 0);

    std::complex<double> intensityTree =
        intensitiesTree->values().at(counter++);

    BOOST_CHECK_EQUAL(intensityNoTree, intensityTree);
    LOG(DEBUG) << "point = " << x << " intensity = " << intensityNoTree
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
  ComPWA::Physics::IntensityBuilderXML Builder;
  auto kin = Builder.createHelicityKinematics(
      partL, tr.get_child("HelicityKinematics"));

  // Generate phsp sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      kin->getParticleStateTransitionKinematicsInfo(), 123));
  std::shared_ptr<ComPWA::Data::DataSet> phspsample(
      ComPWA::Tools::generatePhsp(10000, gen));

  modelStream.clear();
  tr = boost::property_tree::ptree();
  modelStream << HelicityTestModel;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);

  Builder = ComPWA::Physics::IntensityBuilderXML(phspsample);
  auto intens = Builder.createIntensity(partL, kin, tr.get_child("Intensity"));

  ParameterList list;
  intens->addUniqueParametersTo(list);
  LOG(INFO) << "List of parameters: ";
  for (auto p : list.doubleParameters())
    LOG(INFO) << p->to_str();
  BOOST_CHECK_EQUAL(list.doubleParameters().size(), 9);

  // Generate sample
  std::shared_ptr<ComPWA::Data::DataSet> sample(
      ComPWA::Tools::generatePhsp(20, gen));

  sample->convertEventsToParameterList(kin);
  phspsample->convertEventsToParameterList(kin);

  // Testing function tree
  auto tree = intens->createFunctionTree(sample->getParameterList(), "");
  auto tmp = tree->parameter();
  LOG(INFO) << tree->print();
  // Intensity calculated using function tree
  auto intensitiesTree =
      std::dynamic_pointer_cast<Value<std::vector<double>>>(tmp);

  unsigned int counter(0);
  LOG(INFO) << "Loop over events....";
  for (auto const &x : sample->getDataPointList()) {
    // Intensity without function tree
    auto intensityNoTree = intens->evaluate(x);

    double intensityTree = intensitiesTree->values().at(counter++);

    //// Not available in boost 1.54
    // BOOST_TEST(intensityNoTree == intensityTree,
    //           boost::test_tools::tolerance(0.000001));
    BOOST_CHECK_CLOSE(intensityNoTree, intensityTree, 0.0000001);
    LOG(INFO) << "point = " << x << " intensity = " << intensityNoTree
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
  ComPWA::Physics::IntensityBuilderXML Builder;
  auto kin = Builder.createHelicityKinematics(
      partL, tr.get_child("HelicityKinematics"));

  modelStream.clear();
  tr = boost::property_tree::ptree();
  modelStream << SequentialAmplitudeTestModel;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);

  auto seqAmp = Builder.createAmplitude(partL, kin, tr.get_child("Amplitude"));

  ParameterList list;
  seqAmp->addUniqueParametersTo(list);
  LOG(INFO) << "List of parameters: ";
  for (auto p : list.doubleParameters())
    LOG(INFO) << p->to_str();
  BOOST_CHECK_EQUAL(list.doubleParameters().size(), 8);

  // Generate sample
  std::shared_ptr<ComPWA::Generator> gen(new ComPWA::Tools::RootGenerator(
      kin->getParticleStateTransitionKinematicsInfo(), 123));
  std::shared_ptr<ComPWA::Data::DataSet> sample(
      ComPWA::Tools::generatePhsp(20, gen));

  sample->convertEventsToParameterList(kin);

  // Testing function tree
  auto tree = seqAmp->createFunctionTree(sample->getParameterList(), "");
  auto tmp = tree->parameter();
  LOG(INFO) << tree->print();
  // Intensity calculated using function tree
  auto intensitiesTree =
      std::dynamic_pointer_cast<Value<std::vector<std::complex<double>>>>(tmp);

  unsigned int counter(0);
  LOG(INFO) << "Loop over events....";
  for (auto const &x : sample->getDataPointList()) {
    // Intensity without function tree
    auto intensityNoTree = seqAmp->evaluate(x);

    std::complex<double> intensityTree =
        intensitiesTree->values().at(counter++);

    //// Not available in boost 1.54
    // BOOST_TEST(intensityNoTree.real() == intensityTree.real(),
    //           boost::test_tools::tolerance(0.000001));
    // BOOST_TEST(intensityNoTree.imag() == intensityTree.imag(),
    //           boost::test_tools::tolerance(0.000001));
    BOOST_CHECK_CLOSE(intensityNoTree.real(), intensityTree.real(), 0.0000001);
    BOOST_CHECK_CLOSE(intensityNoTree.imag(), intensityTree.imag(), 0.0000001);

    LOG(INFO) << "point = " << x << " intensity = " << intensityNoTree
              << " intensity tree = " << intensityTree;
  }
};

BOOST_AUTO_TEST_SUITE_END()
