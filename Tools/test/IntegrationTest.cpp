// Copyright (c) 2019 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE IntegrationTest

#include "Tools/Integration.hpp"
#include "Core/Logging.hpp"
#include "Data/DataSet.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"
#include "Physics/BuilderXML.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"

#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

using namespace ComPWA;

BOOST_AUTO_TEST_SUITE(ToolTests)

const std::string TestParticles = R"####(
<ParticleList>
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
      <FormFactor Type="BlattWeisskopf"/>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_jpsi'>
        <Value>2.5</Value>
        <Fix>true</Fix>
        <Min>2.0</Min>
        <Max>3.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name='f0'>
    <Pid>6666</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_f0'>
      <Value>1.500</Value>
      <Fix>true</Fix>
      <Min>0.5</Min>
      <Max>1.5</Max>
      <Error>0.00012</Error>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='1'/>
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
</ParticleList>
)####";

const std::string JpsiDecayKinematics = R"####(
<HelicityKinematics>
  <PhspVolume>1</PhspVolume>
  <InitialState>
    <Particle Name='jpsi' PositionIndex='0'/>
  </InitialState>
  <FinalState>
    <Particle Name='gamma' Id='0'/>
    <Particle Name='f0' Id = '1'/>
  </FinalState>
</HelicityKinematics>
)####";

const std::string JpsiDecayTree = R"####(
<Intensity Class='CoherentIntensity'>
  <Amplitude Class='CoefficientAmplitude'>
    <Amplitude Class='SequentialAmplitude'>
      <Amplitude Class='HelicityDecay'>
        <DecayParticle Name='jpsi' Helicity='+1' />
        <DecayProducts>
          <Particle Name='gamma' FinalState='0' Helicity='1' />
          <Particle Name='f0' FinalState='1' Helicity='0' />
        </DecayProducts>
        <CanonicalSum L='0' S='1.0'>
          <ClebschGordan Type='LS' j1='0.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
          <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
        </CanonicalSum>
      </Amplitude>
    </Amplitude>
    <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_0.0_S_1.0;'>
      <Value>1.0</Value>
      <Fix>false</Fix>
    </Parameter>
    <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_0.0_S_1.0;'>
      <Value>0.0</Value>
      <Fix>false</Fix>
    </Parameter>
  </Amplitude>

  <Amplitude Class='CoefficientAmplitude'>
    <Amplitude Class='SequentialAmplitude'>
      <Amplitude Class='HelicityDecay'>
        <DecayParticle Name='jpsi' Helicity='+1' />
        <DecayProducts>
          <Particle Name='gamma' FinalState='0' Helicity='1' />
          <Particle Name='f0' FinalState='1' Helicity='0' />
        </DecayProducts>
        <CanonicalSum L='2' S='1.0'>
          <ClebschGordan Type='LS' j1='2.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
          <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
        </CanonicalSum>
      </Amplitude>
    </Amplitude>
    <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_2.0_S_1.0;'>
      <Value>1.0</Value>
      <Fix>false</Fix>
    </Parameter>
    <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_2.0_S_1.0;'>
      <Value>0.0</Value>
      <Fix>false</Fix>
    </Parameter>
  </Amplitude>
</Intensity>
)####";

// reference gaussian intensity
class Gaussian : public ComPWA::Intensity {
public:
  Gaussian(double Mean, double Width)
      : Mean{"mean", Mean}, Width{"width", Width}, Strength{"strength", 1.0} {}

  std::vector<double> evaluate(const ComPWA::DataMap &Data) noexcept {
    auto const &XValues = Data.at("x");
    std::vector<double> Result(XValues.size());
    std::transform(
        XValues.begin(), XValues.end(), Result.begin(), [&](double x) {
          return 1 / (std::sqrt(2 * M_PI) * Width.Value) * Strength.Value *
                 std::exp(-0.5 * std::pow(x - Mean.Value, 2) /
                          std::pow(Width.Value, 2));
        });
    return Result;
  }

  void updateParametersFrom(const std::vector<double> &Parameters) {
    if (Parameters.size() != 3)
      throw std::runtime_error(
          "IntegrationTest_Gaussian::updateParametersFrom(): Parameter "
          "list size is incorrect!");
    Mean.Value = Parameters[0];
    Width.Value = Parameters[1];
    Strength.Value = Parameters[2];
  }
  std::vector<ComPWA::Parameter> getParameters() const {
    return {Mean, Width, Strength};
  }

private:
  ComPWA::Parameter Mean;
  ComPWA::Parameter Width;
  ComPWA::Parameter Strength;
};

BOOST_AUTO_TEST_CASE(IntegrationAmplitudeModelTest) {
  ComPWA::Logging Log("trace");

  std::stringstream ModelStream;

  // Particle list
  ModelStream << TestParticles;
  auto Particles = readParticles(ModelStream);

  // Kinematics
  ModelStream.clear();
  boost::property_tree::ptree ModelTree;
  ModelStream << JpsiDecayKinematics;
  boost::property_tree::xml_parser::read_xml(ModelStream, ModelTree);
  auto Kinematics = ComPWA::Physics::createHelicityKinematics(
      Particles, ModelTree.get_child("HelicityKinematics"));

  // Estimated using 10^7 events: 0.25672 +- 7.26061e-05
  double TrueIntegral = 0.25672;

  // Generate phsp sample
  ComPWA::Data::Root::RootGenerator Generator(
      Kinematics.getParticleStateTransitionKinematicsInfo());

  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(233);
  auto PhspSample =
      ComPWA::Data::generatePhsp(1000000, Generator, RandomGenerator);

  // Model intensity
  ModelStream.clear();
  ModelTree = boost::property_tree::ptree();
  ModelStream << JpsiDecayTree;
  boost::property_tree::xml_parser::read_xml(ModelStream, ModelTree);

  ComPWA::Physics::IntensityBuilderXML Builder(
      Particles, Kinematics, ModelTree.get_child("Intensity"));
  auto ModelIntensity = Builder.createIntensity();
  auto PhspDataSet = Kinematics.convert(PhspSample);
  auto IntegralSeed1 =
      ComPWA::Tools::integrateWithError(ModelIntensity, PhspDataSet);
  BOOST_CHECK_SMALL(std::abs(IntegralSeed1.first - TrueIntegral),
                    3 * IntegralSeed1.second);

  // Low statistics sample
  ComPWA::Data::resize(PhspDataSet, 200000);
  auto IntegralSeed1Low =
      ComPWA::Tools::integrateWithError(ModelIntensity, PhspDataSet);
  BOOST_CHECK_SMALL(std::abs(IntegralSeed1Low.first - TrueIntegral),
                    3 * IntegralSeed1Low.second);

  // Statistically independent sample
  auto PhspSample2 =
      ComPWA::Data::generatePhsp(1000000, Generator, RandomGenerator);
  auto PhspDataSet2 = Kinematics.convert(PhspSample2);
  auto IntegralSeed2 =
      ComPWA::Tools::integrateWithError(ModelIntensity, PhspDataSet2);
  BOOST_CHECK_SMALL(std::abs(IntegralSeed1.first - IntegralSeed2.first),
                    3 * IntegralSeed1.second);
  BOOST_CHECK_SMALL(std::abs(IntegralSeed2.first - TrueIntegral),
                    3 * IntegralSeed2.second);
}

BOOST_AUTO_TEST_CASE(IntegrationGaussianTest) {
  ComPWA::Logging Log("trace");

  // the reference normal distribution
  double Mean(3.0);
  double Sigma(0.1);

  auto Gauss = Gaussian(Mean, Sigma);

  std::pair<double, double> DomainRange(Mean - 10.0 * Sigma,
                                        Mean + 10.0 * Sigma);
  ComPWA::Data::DataSet PhspSample;
  PhspSample.Data.insert(std::make_pair("x", std::vector<double>()));
  std::mt19937 MtGen(123456);

  std::uniform_real_distribution<double> Distribution(DomainRange.first,
                                                      DomainRange.second);

  for (unsigned int i = 0; i < 200000; ++i) {
    PhspSample.Data["x"].push_back(Distribution(MtGen));
    PhspSample.Weights.push_back(1.0);
  }

  auto Integral = ComPWA::Tools::integrateWithError(
      Gauss, PhspSample, DomainRange.second - DomainRange.first);
  LOG(INFO) << "Calculated integral: " << Integral.first << "+-"
            << Integral.second;

  BOOST_CHECK_SMALL(std::abs(Integral.first - 1.0), 3 * Integral.second);
}

BOOST_AUTO_TEST_SUITE_END()
