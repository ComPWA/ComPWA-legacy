// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE FitTest

#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"
#include "Physics/BuilderXML.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Plotting/RootPlotData.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>

using namespace ComPWA;
using ComPWA::Optimizer::Minuit2::MinuitResult;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

// We define an intensity model using a raw string literal. Currently, this is
// just a toy model without any physical meaning.
// (comments within the string are ignored!). This is convenient since we
// do not have to configure the build system to copy input files somewhere.
// In practice you may want to use a normal XML input file instead.
std::string AmplitudeModel = R"####(
<Intensity Class='NormalizedIntensity'>
  <Intensity Class='CoherentIntensity'>
    <Amplitude Class="CoefficientAmplitude">
      <Parameter Class='Double' Type="Magnitude"  Name="Magnitude_f2">
        <Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>5.0</Max>
        <Fix>false</Fix>
      </Parameter>
      <Parameter Class='Double' Type="Phase" Name="Phase_f2">
        <Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
        <Fix>false</Fix>
      </Parameter>
	  <Amplitude Class="SequentialAmplitude">
      <Amplitude Class="HelicityDecay">
        <DecayParticle Name="J/psi" Helicity="0"/>
        <DecayProducts>
          <Particle Name="f2(1270)" FinalState="1 2"  Helicity="0"/>
          <Particle Name="gamma" FinalState="0"  Helicity="1"/>
        </DecayProducts>
      </Amplitude>
      <Amplitude Class="HelicityDecay">
        <DecayParticle Name="f2(1270)" Helicity="0"/>
        <RecoilSystem FinalState="0" />
        <DecayProducts>
          <Particle Name="pi0" FinalState="1"  Helicity="0"/>
          <Particle Name="pi0" FinalState="2"  Helicity="0"/>
        </DecayProducts>
      </Amplitude>
      </Amplitude>
    </Amplitude>
    <Amplitude Class="CoefficientAmplitude">
      <Parameter Class='Double' Type="Magnitude"  Name="Magnitude_my">
        <Value>3.0</Value>
        <Min>-1.0</Min>
        <Max>5.0</Max>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type="Phase" Name="Phase_my`">
        <Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
        <Fix>true</Fix>
      </Parameter>
      <Amplitude Class="SequentialAmplitude">
      <Amplitude Class="HelicityDecay">
        <DecayParticle Name="J/psi" Helicity="0"/>
        <DecayProducts>
          <Particle Name="myRes" FinalState="1 2"  Helicity="0"/>
          <Particle Name="gamma" FinalState="0"  Helicity="1"/>
        </DecayProducts>
      </Amplitude>
      <Amplitude Class="HelicityDecay">
        <DecayParticle Name="myRes" Helicity="0"/>
        <RecoilSystem FinalState="0" />
        <DecayProducts>
          <Particle Name="pi0" FinalState="1"  Helicity="0"/>
          <Particle Name="pi0" FinalState="2"  Helicity="0"/>
        </DecayProducts>
      </Amplitude>
      </Amplitude>
    </Amplitude>
  </Intensity>
</Intensity>
)####";

std::string MyParticleList = R"####(
<ParticleList>
	<Particle Name="J/psi">
		<Pid>443</Pid>
		<Parameter Class="Double" Type="Mass" Name="Mass_jpsi">
			<Value>3.096900</Value>
			<Fix>true</Fix>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="1"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="-1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="-1"/>
		<DecayInfo Type="nonResonant">
		</DecayInfo>
	</Particle>
  <Particle Name="f2(1270)">
    <Pid>225</Pid>
    <Parameter Class='Double' Type="Mass" Name="Mass_f2(1270)">
      <Value>1.2755</Value>
      <Error>8.0E-04</Error>
      <Min>0.1</Min>
      <Max>2.0</Max>
      <Fix>false</Fix>
    </Parameter>
    <QuantumNumber Class="Spin" Type="Spin" Value="2"/>
    <QuantumNumber Class="Int" Type="Charge" Value="0"/>
    <QuantumNumber Class="Int" Type="Parity" Value="+1"/>
    <QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
    <DecayInfo Type="relativisticBreitWignerAC">
      <Parameter Class='Double' Type="Width" Name="Width_f2(1270)">
        <Value>0.1867</Value>
      </Parameter>
      <Parameter Class='Double' Type="MesonRadius" Name="Radius_rho">
        <Value>2.5</Value>
        <Fix>true</Fix>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name="myRes">
    <Pid>999999</Pid>
    <Parameter Class='Double' Type="Mass" Name="Mass_myRes">
      <Value>2.0</Value>
      <Error>8.0E-04</Error>
      <Min>1.1</Min>
      <Max>4.0</Max>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class="Spin" Type="Spin" Value="1"/>
    <QuantumNumber Class="Int" Type="Charge" Value="0"/>
    <QuantumNumber Class="Int" Type="Parity" Value="+1"/>
    <QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
    <DecayInfo Type="relativisticBreitWignerAC">
      <Parameter Class='Double' Type="Width" Name="Width_myRes">
        <Value>1.0</Value>
        <Min>0.1</Min>
        <Max>1.0</Max>
        <Fix>false</Fix>
      </Parameter>
      <Parameter Class='Double' Type="MesonRadius" Name="Radius_myRes">
        <Value>2.5</Value>
        <Fix>true</Fix>
      </Parameter>
    </DecayInfo>
  </Particle>
</ParticleList>
)####";

FitParameter<double> getFitParameter(FitParameterList ParameterList,
                                     std::string Name) {
  auto Result =
      std::find_if(ParameterList.begin(), ParameterList.end(),
                   [&Name](const ComPWA::FitParameter<double> &Parameter) {
                     return Parameter.Name == Name;
                   });
  if (Result == ParameterList.end())
    return FitParameter<double>();
  return *Result;
}

BOOST_AUTO_TEST_SUITE(FitTest)

BOOST_AUTO_TEST_CASE(HelicityDalitzFit) {
  ComPWA::Logging Log("debug");

  std::stringstream ParticlesStream;
  ParticlesStream << MyParticleList;
  // List with all particle information needed
  auto ParticleList = readParticles("particle_list.xml");
  insertParticles(ParticleList, ParticlesStream);

  //---------------------------------------------------
  // 1) Create Kinematics object
  //---------------------------------------------------
  std::vector<pid> InitialState = {443};
  std::vector<pid> FinalState = {22, 111, 111};
  HelicityKinematics Kinematics(ParticleList, InitialState, FinalState);

  //---------------------------------------------------
  // 2) Generate a large phase space sample
  //---------------------------------------------------
  ComPWA::Data::Root::RootGenerator Generator(
      Kinematics.getParticleStateTransitionKinematicsInfo());

  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(173);

  auto PhspSample =
      ComPWA::Data::generatePhsp(100000, Generator, RandomGenerator);

  //---------------------------------------------------
  // 3) Create intensity from pre-defined model
  //---------------------------------------------------
  // Read in model property_tree
  std::stringstream ModelStream;
  ModelStream << AmplitudeModel;
  boost::property_tree::ptree ModelTree;
  boost::property_tree::xml_parser::read_xml(ModelStream, ModelTree);

  // Construct intensity class from model string
  ComPWA::Physics::IntensityBuilderXML Builder(
      ParticleList, Kinematics, ModelTree.get_child("Intensity"), PhspSample);
  auto Intensity = Builder.createIntensity();

  //---------------------------------------------------
  // 4) Generate a data sample given intensity and kinematics
  //---------------------------------------------------
  RandomGenerator.setSeed(1234);

  auto DataSample = ComPWA::Data::generate(1000, Kinematics, RandomGenerator,
                                           Intensity, PhspSample);

  auto PhspSampleDataSet = Kinematics.convert(PhspSample);
  auto SampleDataSet = Kinematics.convert(DataSample);

  //---------------------------------------------------
  // 5) Fit the model to the data and print the result
  //---------------------------------------------------
  auto Estimator = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
      Intensity, SampleDataSet);

  auto Optimizer = Optimizer::Minuit2::MinuitIF();

  auto FitParams = std::get<1>(Estimator);

  // STARTING MINIMIZATION
  ComPWA::Optimizer::Minuit2::MinuitResult FitResult =
      Optimizer.optimize(std::get<0>(Estimator), FitParams);

  LOG(INFO) << FitResult;

  // output << result->finalLH();
  BOOST_CHECK_EQUAL(DataSample.Events.size(), 1000);
  BOOST_CHECK_CLOSE(FitResult.FinalEstimatorValue, -730, 5.); // 5% tolerance
  double Sigma(3.0);

  auto FitParameter =
      getFitParameter(FitResult.FinalParameters, "Magnitude_f2");
  BOOST_CHECK_GT(FitParameter.Value + Sigma * FitParameter.Error.second, 0.0);
  BOOST_CHECK_GT(1.0, FitParameter.Value - Sigma * FitParameter.Error.first);

  FitParameter = getFitParameter(FitResult.FinalParameters, "Phase_f2");
  BOOST_CHECK_GT((FitParameter.Value + Sigma * FitParameter.Error.second), 0.0);
  BOOST_CHECK_GT(0.0, (FitParameter.Value - Sigma * FitParameter.Error.first));
  FitParameter = getFitParameter(FitResult.FinalParameters, "Mass_f2(1270)");
  BOOST_CHECK_GT((FitParameter.Value + Sigma * FitParameter.Error.second),
                 1.2755);
  BOOST_CHECK_GT(1.2755,
                 (FitParameter.Value - Sigma * FitParameter.Error.first));
  FitParameter = getFitParameter(FitResult.FinalParameters, "Width_myRes");
  BOOST_CHECK_GT((FitParameter.Value + Sigma * FitParameter.Error.second), 1.0);
  BOOST_CHECK_GT(1.0, (FitParameter.Value - Sigma * FitParameter.Error.first));
};
BOOST_AUTO_TEST_SUITE_END()
