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
#include "Core/FitParameter.hpp"
#include "Data/DataSet.hpp"
#include "Physics/Amplitude.hpp"
#include "Physics/Dynamics/RelativisticBreitWigner.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Physics/Dynamics/FormFactorDecorator.hpp"
#include "Physics/Dynamics/NonResonant.hpp"
#include "Physics/Dynamics/RelativisticBreitWigner.hpp"
#include "Tools/Generate.hpp"
#include "Tools/RootGenerator.hpp"

#include <boost/foreach.hpp>
#include <boost/locale/utf.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

using namespace ComPWA::Physics::Dynamics;
using namespace ComPWA::Physics;

BOOST_AUTO_TEST_SUITE(Dynamics)

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
    <DecayInfo Type='nonResonant'>
      <FormFactor Type='1' />
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
      <FormFactor Type='1' />
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

const std::string JpsiDecayKinematics = R"####(
<HelicityKinematics>
  <PhspVolume>1</PhspVolume>
  <InitialState>
    <Particle Name='jpsi' PositionIndex='0'/>
  </InitialState>
  <FinalState>
    <Particle Name='pi0' Id='0'/>
    <Particle Name='omega' Id='1'/>
  </FinalState>
</HelicityKinematics>
)####";
//jpsi -> omega pi0
//         L 0 S h   J h  s2 h2 s3 -h3 S h
//         1/sqrt{2}        1
//F_0_1 = (1 0 1 -1|1 -1) (0 0 1 -1|1 -1) a_1_1
const std::string JpsiDecayTree = R"####(
<Amplitude Class='HelicityDecay' Name='jpsitoPi0Omega'>
  <DecayParticle Name='jpsi' Helicity='+1' />
  <DecayProducts>
    <Particle Name='pi0' FinalState='0' Helicity='0' />
    <Particle Name='omega' FinalState='1' Helicity='+1' />
  </DecayProducts>
  <CanonicalSum L='1.0' S='1.0'>
    <ClebschGordan Type='LS' j1='1.0' m1='0.0' j2='1.0' m2='-1.0' J='1.0' M='-1.0' />
    <ClebschGordan Type='s2s3' j1='0.0' m1='0.0' j2='1.0' m2='-1.0' J='1.0' M='-1.0' />
  </CanonicalSum>
</Amplitude>
)####";

const std::string OmegaDecayKinematics = R"####(
<HelicityKinematics>
  <PhspVolume>1</PhspVolume>
  <InitialState>
    <Particle Name='omega' PositionIndex='0'/>
  </InitialState>
  <FinalState>
    <Particle Name='gamma' Id='0'/>
    <Particle Name='pi0' Id='1'/>
  </FinalState>
</HelicityKinematics>
)####";
const std::string OmegaDecayTree = R"####(
<Amplitude Class='HelicityDecay' Name='omegatoPi0Gamma'>
  <DecayParticle Name='omega' Helicity='+1' />
  <DecayProducts>
    <Particle Name='pi0' FinalState='0' Helicity='0' />
    <Particle Name='gamma' FinalState='1' Helicity='+1' />
  </DecayProducts>
  <CanonicalSum L='1.0' S='1.0'>
    <ClebschGordan Type='LS' j1='1.0' m1='0.0' j2='1.0' m2='-1.0' J='1.0' M='-1.0'/>
    <ClebschGordan Type='s2s3' j1='0.0' m1='0.0' j2='1.0' m2='-1.0' J='1.0' M='-1.0'/>
  </CanonicalSum>
</Amplitude>
)####";

void testFormFactorDecoratedBW(const std::string &partList,
    const std::string &kinematics,
    const std::string &decayTree, const std::string &mother,
    const std::string &daughter1, const std::string &daughter2);

// IF (BW->evaluate() * FormFactor) == (FormFactorDecorator(BW)->evaluate())
//jpsi -> omega pi0 and omega -> gamma pi0
//         LS              s2s3
//F_0_1 = (1 0 1 -1|1 -1) (0 0 1 -1|1 -1) a_1_1
BOOST_AUTO_TEST_CASE(BWWithFormFactorDecorator) {
  ComPWA::Logging log("", "trace");

  LOG(INFO) << "Now check formFactorDecorator for nonResonant...";
  //test Undecorated NonResonant times Production Formfactor
  //and Production Formfactor decorated NonResonant
  testFormFactorDecoratedBW(HelicityTestParticles, JpsiDecayKinematics,
      JpsiDecayTree, "jpsi", "omega", "pi0");

  //test Undecorated RelativisticBreitWigner times Production FormFactor
  //and Production Formfactor decorated RelativisticBreitWigner
  LOG(INFO) << "Now check formFactorDecorator for relativisticBreitWigner...";
  testFormFactorDecoratedBW(HelicityTestParticles, OmegaDecayKinematics,
      OmegaDecayTree, "omega", "gamma", "pi0");

}

void testFormFactorDecoratedBW(const std::string &partList,
    const std::string &kinematics,
    const std::string &decayTree, const std::string &mother,
    const std::string &daughter1, const std::string &daughter2) {
  boost::property_tree::ptree tr;
  std::stringstream modelStream;
  // Construct particle list from XML tree
  modelStream << partList;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, tr);

  ComPWA::Physics::IntensityBuilderXML Builder;

  //check nonResonant with jpsi -> omega pi0
  modelStream.clear();
  tr = boost::property_tree::ptree();
  modelStream << kinematics;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);

  auto jpsiKin = Builder.createHelicityKinematics(partL,
      tr.get_child("HelicityKinematics"));

  // Generate sample                                                                
  auto jpsiGen = std::make_shared<ComPWA::Tools::RootGenerator>(
      jpsiKin->getParticleStateTransitionKinematicsInfo(), 123);
  auto jpsiSample = ComPWA::Tools::generatePhsp(20, jpsiGen);

  // add subsystem to kinematics
  modelStream.clear();
  tr = boost::property_tree::ptree();
  modelStream.clear();
  modelStream << JpsiDecayTree;
  boost::property_tree::xml_parser::read_xml(modelStream, tr);
  unsigned int subSysIndex(jpsiKin->addSubSystem(
      SubSystem(tr.get_child("Amplitude"))));
  unsigned int dataPos = 3 * subSysIndex;
 
  jpsiSample->convertEventsToDataPoints(jpsiKin);

  const std::string BWTag = partL->find(mother)->second.GetDecayInfo()
      .get<std::string>("<xmlattr>.Type");

  std::shared_ptr<ComPWA::Physics::Dynamics::AbstractDynamicalFunction>
      jpsiNoProdFF;

  std::string decayName = mother + "_to_" + daughter1 + "+" + daughter2;
  //Undecorated NonResonant
  if (BWTag == "nonResonant") {
    jpsiNoProdFF = std::make_shared<ComPWA::Physics::Dynamics::NonResonant>(
        mother);
  } else if (BWTag == "relativisticBreitWigner") {
    jpsiNoProdFF = std::make_shared<ComPWA::Physics::Dynamics
        ::RelativisticBreitWigner>(mother,
        std::pair<std::string, std::string>(daughter1, daughter2),
        partL);
  }

  // set formfactor decorator
  int ffType = 1;
  int orbitL = 1;
  auto parMass1 = std::make_shared<ComPWA::FitParameter>(                                
      partL->find(daughter1)->second.GetMassPar());                              
  auto parMass2 = std::make_shared<ComPWA::FitParameter>(                                
      partL->find(daughter2)->second.GetMassPar());                              
  std::shared_ptr<ComPWA::FitParameter> parRadius;
  const auto &decayInfo = partL->find(mother)->second.GetDecayInfo();
          
  for (auto const &node : decayInfo.get_child("")) {
    if (node.first != "Parameter") continue;
    if (node.second.get<std::string>("<xmlattr>.Type") != "MesonRadius")
      continue;
    parRadius = std::make_shared<ComPWA::FitParameter>(node.second);
  }

  FormFactorDecorator jpsiProdFF(mother,
      jpsiNoProdFF, parMass1, parMass2, parRadius,
      (ComPWA::Spin) orbitL, (FormFactorType) ffType);
  
  //compare FormFactorDecoratedNonResonant 
  //and UndecoratedNonResonant * FormFactor
  for (const auto &point : jpsiSample->getDataPointList()) {
    std::complex<double> value = jpsiNoProdFF->evaluate(point, dataPos);
    //calculate production formfactor by hand
    double ff = FormFactor(sqrt(point.KinematicVariableList[dataPos]),
        parMass1->value(), parMass2->value(), (unsigned int) orbitL,
        parRadius->value(), (FormFactorType) ffType);
    std::complex<double> valueTimesProdFF = value * ff;

    //formfactor decorated BW
    std::complex<double> valueProdFFDecorated = 
        jpsiProdFF.evaluate(point, dataPos);

    LOG(INFO) << "<  BreitWignerType: " << BWTag << " L = " << orbitL
        << " FormFactorType = " << ffType;
    LOG(INFO) << "     BWFormFactorDecorated = " << valueProdFFDecorated << ",";
    LOG(INFO) << "FormFactor * BWUndecorated = " << valueTimesProdFF << " >";
    BOOST_CHECK_EQUAL(valueProdFFDecorated, valueTimesProdFF);
  }

  // compare Tree->Parameter and Evaluate
  jpsiSample->convertEventsToParameterList(jpsiKin);
  auto funcTree = jpsiProdFF.createFunctionTree(jpsiSample->getParameterList(),
      dataPos, "");
  auto tmp = funcTree->parameter();
  LOG(INFO) << funcTree->print();
  auto treeValues =                                                         
      std::dynamic_pointer_cast<ComPWA::Value<std::vector
      <std::complex<double>>>>(tmp); 
  unsigned int pointIndex = 0;
  for (const auto &point : jpsiSample->getDataPointList()) {
    std::complex<double> valueProdFFDecorated = 
        jpsiProdFF.evaluate(point, dataPos);
    std::complex<double> valueFromTreeParameter
        = treeValues->values().at(pointIndex++);
    LOG(INFO) << "<  BreitWignerType: " << BWTag << " L = " << orbitL
        << " FormFactorType = " << ffType;
    LOG(INFO) << "       BWFormFactorDecorated.evaluate() = " 
        << valueProdFFDecorated << ",";
    LOG(INFO) << " BWFormFactorDecorated.Tree.Parameter() = " 
        << valueFromTreeParameter << " >";
    BOOST_CHECK_EQUAL(valueProdFFDecorated, valueFromTreeParameter);
  }
}


BOOST_AUTO_TEST_SUITE_END()
