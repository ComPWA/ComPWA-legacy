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
std::string amplitudeModel = R"####(
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

std::string myParticles = R"####(
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

FitParameter<double> getFitParameter(FitParameterList list, std::string name) {
  auto res = std::find_if(list.begin(), list.end(),
                          [&name](const ComPWA::FitParameter<double> &x) {
                            return x.Name == name;
                          });
  if (res == list.end())
    return FitParameter<double>();
  return *res;
}

BOOST_AUTO_TEST_SUITE(FitTest)

BOOST_AUTO_TEST_CASE(HelicityDalitzFit) {
  ComPWA::Logging log("debug");

  std::stringstream ParticlesStream;
  ParticlesStream << myParticles;
  // List with all particle information needed
  auto partL = readParticles("particle_list.xml");
  insertParticles(partL, ParticlesStream);

  //---------------------------------------------------
  // 1) Create Kinematics object
  //---------------------------------------------------
  std::vector<pid> initialState = {443};
  std::vector<pid> finalState = {22, 111, 111};
  HelicityKinematics kin(partL, initialState, finalState);

  //---------------------------------------------------
  // 2) Generate a large phase space sample
  //---------------------------------------------------
  ComPWA::Data::Root::RootGenerator gen(
      kin.getParticleStateTransitionKinematicsInfo());

  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(173);

  auto phspSample = ComPWA::Data::generatePhsp(100000, gen, RandomGenerator);

  //---------------------------------------------------
  // 3) Create intensity from pre-defined model
  //---------------------------------------------------
  // Read in model property_tree
  std::stringstream modelStream;
  modelStream << amplitudeModel;
  boost::property_tree::ptree modelTree;
  boost::property_tree::xml_parser::read_xml(modelStream, modelTree);

  // Construct intensity class from model string
  ComPWA::Physics::IntensityBuilderXML Builder(
      partL, kin, modelTree.get_child("Intensity"), phspSample);
  auto intens = Builder.createIntensity();

  //---------------------------------------------------
  // 4) Generate a data sample given intensity and kinematics
  //---------------------------------------------------
  RandomGenerator.setSeed(1234);

  auto sample =
      ComPWA::Data::generate(1000, kin, RandomGenerator, intens, phspSample);

  auto PhspSampleDataSet = kin.convert(phspSample);
  auto SampleDataSet = kin.convert(sample);

  //---------------------------------------------------
  // 5) Fit the model to the data and print the result
  //---------------------------------------------------
  auto esti = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
      intens, SampleDataSet);

  // LOG(INFO) << esti.first.print(25);

  auto minuitif = Optimizer::Minuit2::MinuitIF();

  auto FitParams = std::get<1>(esti);

  // STARTING MINIMIZATION
  ComPWA::Optimizer::Minuit2::MinuitResult result =
      minuitif.optimize(std::get<0>(esti), FitParams);

  LOG(INFO) << result;

  // output << result->finalLH();
  BOOST_CHECK_EQUAL(sample.size(), 1000);
  BOOST_CHECK_CLOSE(result.FinalEstimatorValue, -730, 5.); // 5% tolerance
  double sigma(3.0);

  auto fitpar = getFitParameter(result.FinalParameters, "Magnitude_f2");
  BOOST_CHECK_GT(fitpar.Value + sigma * fitpar.Error.second, 0.0);
  BOOST_CHECK_GT(1.0, fitpar.Value - sigma * fitpar.Error.first);

  fitpar = getFitParameter(result.FinalParameters, "Phase_f2");
  BOOST_CHECK_GT((fitpar.Value + sigma * fitpar.Error.second), 0.0);
  BOOST_CHECK_GT(0.0, (fitpar.Value - sigma * fitpar.Error.first));
  fitpar = getFitParameter(result.FinalParameters, "Mass_f2(1270)");
  BOOST_CHECK_GT((fitpar.Value + sigma * fitpar.Error.second), 1.2755);
  BOOST_CHECK_GT(1.2755, (fitpar.Value - sigma * fitpar.Error.first));
  fitpar = getFitParameter(result.FinalParameters, "Width_myRes");
  BOOST_CHECK_GT((fitpar.Value + sigma * fitpar.Error.second), 1.0);
  BOOST_CHECK_GT(1.0, (fitpar.Value - sigma * fitpar.Error.first));
};
BOOST_AUTO_TEST_SUITE_END()
