// Copyright (c) 2019 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE IntegrationTest

#include "Tools/Integration.hpp"
#include "Core/Logging.hpp"
#include "Data/DataSet.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"

#include <boost/test/unit_test.hpp>

using namespace ComPWA;

BOOST_AUTO_TEST_SUITE(ToolTests)

const std::string TestParticles = R"####(
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
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='1' />
      <Parameter Class='Double' Type='Width' Name='Width_f0'>
        <Value>0.050</Value>
        <Fix>true</Fix>
        <Min>0.0</Min>
        <Max>1.0</Max>
        <Error>0.00008</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_f0'>
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
    <Particle Name='gamma' Id='0'/>
    <Particle Name='f0' Id = '1'/>
  </FinalState>
</HelicityKinematics>
)####";

const std::string JpsiDecayTree = R"####(
    <Intensity Class='CoherentIntensity' Name='coherent_0'>
      <Amplitude Class='CoefficientAmplitude' Name='jpsi_1_to_gamma_1+f0_0_L_0.0_S_1.0;'>
        <Amplitude Class='SequentialAmplitude' Name='jpsi_1_to_gamma_1+f0_0_L_0.0_S_1.0;'>
          <Amplitude Class='HelicityDecay' Name='jpsi_to_gamma_f0_L0_S1'>
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

      <Amplitude Class='CoefficientAmplitude' Name='jpsi_1_to_gamma_1+f0_0_L_2.0_S_1.0;'>
        <Amplitude Class='SequentialAmplitude' Name='jpsi_1_to_gamma_1+f0_0_L_2.0_S_1.0;'>
          <Amplitude Class='HelicityDecay' Name='jpsi_to_gamma_f0_L2_S1'>
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
  Gaussian(double mean, double width)
      : Mean{"mean", mean}, Width{"width", width}, Strength{"strength", 1.0} {}

  std::vector<double>
  evaluate(const std::vector<std::vector<double>> &data) noexcept {
    auto const &xvals = data[0];
    std::vector<double> result(xvals.size());
    std::transform(xvals.begin(), xvals.end(), result.begin(), [&](double x) {
      return 1 / (std::sqrt(2 * M_PI) * Width.Value) * Strength.Value *
             std::exp(-0.5 * std::pow(x - Mean.Value, 2) /
                      std::pow(Width.Value, 2));
    });
    return result;
  }

  void updateParametersFrom(const std::vector<double> &params) {
    if (params.size() != 3)
      throw std::runtime_error(
          "IntegrationTest_Gaussian::updateParametersFrom(): Parameter "
          "list size is incorrect!");
    Mean.Value = params[0];
    Width.Value = params[1];
    Strength.Value = params[2];
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
  ComPWA::Logging Log("trace", "");

  boost::property_tree::ptree ModelTree;
  std::stringstream ModelStream;
  ComPWA::Physics::IntensityBuilderXML Builder;

  // Particle list
  ModelStream << TestParticles;
  boost::property_tree::xml_parser::read_xml(ModelStream, ModelTree);
  auto PartL = std::make_shared<ComPWA::PartList>();
  ReadParticles(PartL, ModelTree);

  // Kinematics
  ModelStream.clear();
  ModelTree = boost::property_tree::ptree();
  ModelStream << JpsiDecayKinematics;
  boost::property_tree::xml_parser::read_xml(ModelStream, ModelTree);
  auto DecayKin = Builder.createHelicityKinematics(
      PartL, ModelTree.get_child("HelicityKinematics"));

  // Estimated using 10^7 events: 0.198984549+-5.627720404e-05
  double TrueIntegral = 0.198984549;

  // Generate phsp sample
  ComPWA::Data::Root::RootGenerator Gen(
      DecayKin.getParticleStateTransitionKinematicsInfo());

  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(233);
  auto PhspSample = Data::generatePhsp(1000000, Gen, RandomGenerator);

  Builder = ComPWA::Physics::IntensityBuilderXML(PhspSample);

  // Model intensity
  ModelStream.clear();
  ModelTree = boost::property_tree::ptree();
  ModelStream << JpsiDecayTree;
  boost::property_tree::xml_parser::read_xml(ModelStream, ModelTree);

  auto ModelIntensity = Builder.createIntensity(
      PartL, DecayKin, ModelTree.get_child("Intensity"));

  auto PhspDataSet = ComPWA::Data::convertEventsToDataSet(PhspSample, DecayKin);
  auto integ_seed1 =
      ComPWA::Tools::integrateWithError(ModelIntensity, PhspDataSet);
  BOOST_CHECK_SMALL(std::abs(integ_seed1.first - TrueIntegral),
                    3 * integ_seed1.second);

  // Low statistics sample
  ComPWA::Data::resize(PhspDataSet, 200000);
  auto integ_seed1_low =
      ComPWA::Tools::integrateWithError(ModelIntensity, PhspDataSet);
  BOOST_CHECK_SMALL(std::abs(integ_seed1_low.first - TrueIntegral),
                    3 * integ_seed1_low.second);

  // Statistically independent sample
  auto PhspSample2 = Data::generatePhsp(1000000, Gen, RandomGenerator);
  auto PhspDataSet2 =
      ComPWA::Data::convertEventsToDataSet(PhspSample2, DecayKin);
  auto integ_seed2 =
      ComPWA::Tools::integrateWithError(ModelIntensity, PhspDataSet2);
  BOOST_CHECK_SMALL(std::abs(integ_seed1.first - integ_seed2.first),
                    3 * integ_seed1.second);
  BOOST_CHECK_SMALL(std::abs(integ_seed2.first - TrueIntegral),
                    3 * integ_seed2.second);
}

BOOST_AUTO_TEST_CASE(IntegrationGaussianTest) {
  ComPWA::Logging Log("trace", "");

  // the reference normal distribution
  double mean(3.0);
  double sigma(0.1);

  auto Gauss = Gaussian(mean, sigma);

  std::pair<double, double> domain_range(mean - 10.0 * sigma,
                                         mean + 10.0 * sigma);
  ComPWA::Data::DataSet PhspSample;
  PhspSample.Data.push_back({});
  std::mt19937 mt_gen(123456);

  std::uniform_real_distribution<double> distribution(domain_range.first,
                                                      domain_range.second);

  for (unsigned int i = 0; i < 200000; ++i) {
    PhspSample.Data[0].push_back(distribution(mt_gen));
    PhspSample.Weights.push_back(1.0);
  }

  auto integral = ComPWA::Tools::integrateWithError(
      Gauss, PhspSample, domain_range.second - domain_range.first);
  LOG(INFO) << "Calculated integral: " << integral.first << "+-"
            << integral.second;

  BOOST_CHECK_SMALL(std::abs(integral.first - 1.0), 3 * integral.second);
}

BOOST_AUTO_TEST_SUITE_END()
