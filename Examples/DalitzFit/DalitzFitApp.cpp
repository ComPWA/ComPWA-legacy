// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Simple Dalitz plot analysis with ComPWA
///

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Physics/ParticleList.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Plotting/DalitzPlot.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

using namespace ComPWA;
using ComPWA::Optimizer::Minuit2::MinuitResult;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

// Enable serialization of MinuitResult. For some reason has to be outside
// any namespaces.
// BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)

// We define an intensity model using a raw string literal. Currently, this is
// just a toy model without any physical meaning.
// (comments within the string are ignored!). This is convenient since we
// do not have to configure the build system to copy input files somewhere.
// In practise you may want to use a normal XML input file instead.
std::string amplitudeModel = R"####(
<Intensity Class="NormalizedIntensity" Name="jpsiGammaPiPi_norm">
  <IntegrationStrategy Class="MCIntegrationStrategy"/>
  <Intensity Class="CoherentIntensity" Name="jpsiGammaPiPi">
    <Amplitude Class="CoefficientAmplitude" Name="f2(1270)">
	  <Parameter Class='Double' Type="Magnitude"  Name="Magnitude_f2">
		<Value>1.0</Value>
		<Min>-1.0</Min>
		<Max>2.0</Max>
		<Fix>false</Fix>
	  </Parameter>
	  <Parameter Class='Double' Type="Phase" Name="Phase_f2">
		<Value>0.0</Value>
		<Min>-100</Min>
		<Max>100</Max>
		<Fix>false</Fix>
	  </Parameter>
      <Amplitude Class="NormalizedAmplitude" Name="f2(1270)_normed">
        <IntegrationStrategy Class="MCIntegrationStrategy"/>
		<Amplitude Class="HelicityDecay" Name="f2ToPiPi">
		  <DecayParticle Name="f2(1270)" Helicity="0"/>
		  <RecoilSystem FinalState="0" />
		  <DecayProducts>
			<Particle Name="pi0" FinalState="1"  Helicity="0"/>
			<Particle Name="pi0" FinalState="2"  Helicity="0"/>
		  </DecayProducts>
		</Amplitude>
      </Amplitude>
    </Amplitude>
	<Amplitude Class="CoefficientAmplitude" Name="myAmp">
	  <Parameter Class='Double' Type="Magnitude"  Name="Magnitude_my">
		<Value>1.0</Value>
		<Min>-1.0</Min>
		<Max>2.0</Max>
		<Fix>true</Fix>
	  </Parameter>
	  <Parameter Class='Double' Type="Phase" Name="Phase_my`">
		<Value>0.0</Value>
		<Min>-100</Min>
		<Max>100</Max>
		<Fix>true</Fix>
	  </Parameter>
      <Amplitude Class="NormalizedAmplitude" Name="myAmp_normed">
        <IntegrationStrategy Class="MCIntegrationStrategy"/>
		<Amplitude Class="HelicityDecay" Name="MyResToPiPi">
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
    <DecayInfo Type="relativisticBreitWigner">
      <FormFactor Type="0" />
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
    <DecayInfo Type="relativisticBreitWigner">
      <FormFactor Type="0" />
      <Parameter Class='Double' Type="Width" Name="Width_myRes">
        <Value>1.0</Value>
        <Min>0.1</Min>
        <Max>1.5</Max>
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

///
/// Simple Dalitz plot fit of the channel J/psi -> gamma pi0 pi0
///
/// The basic procedure is the following:
/// 1) Create Kinematics object
/// 2) Generate a large phase space sample
/// 3) Create intensity from pre-defined model
/// 4) Generate a data sample given intensity and kinematics
/// 5) Fit the model to the data
/// 6) Plot data sample and intensity
///
int main(int argc, char **argv) {
  // initialize logging
  Logging log("debug", "DalitzFit-log.txt");

  // List with all particle information needed
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, ComPWA::Physics::defaultParticleList);
  ReadParticles(partL, myParticles);

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

  auto phspSample(ComPWA::Data::generatePhsp(1000000, gen, RandomGenerator));

  //---------------------------------------------------
  // 3) Create intensity from pre-defined model
  //---------------------------------------------------
  // Read in model property_tree
  std::stringstream modelStream;
  modelStream << amplitudeModel;
  boost::property_tree::ptree modelTree;
  boost::property_tree::xml_parser::read_xml(modelStream, modelTree);

  // Construct intensity class from model string
  ComPWA::Physics::IntensityBuilderXML Builder(phspSample);
  auto intens =
      Builder.createIntensity(partL, kin, modelTree.get_child("Intensity"));

  //---------------------------------------------------
  // 4) Generate a data sample given intensity and kinematics
  //---------------------------------------------------

  auto sample =
      ComPWA::Data::generate(1000, kin, RandomGenerator, intens, phspSample);

  //---------------------------------------------------
  // 5) Fit the model to the data and print the result
  //---------------------------------------------------

  auto SampleDataSet = Data::convertEventsToDataSet(sample, kin);

  auto esti = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
      intens, SampleDataSet);

  LOG(DEBUG) << esti.first.print(25);

  auto minuitif = Optimizer::Minuit2::MinuitIF();

  // STARTING MINIMIZATION
  auto result = minuitif.optimize(std::get<0>(esti), std::get<1>(esti));

  // TODO: do fit fraction calculation part

  LOG(INFO) << result;

  //---------------------------------------------------
  // 5.1) Save the fit result
  //---------------------------------------------------
  std::ofstream ofs("DalitzFit-fitResult.xml");
  boost::archive::xml_oarchive oa(ofs);
  oa << BOOST_SERIALIZATION_NVP(result);

  //---------------------------------------------------
  // 6) Plot data sample and intensity
  //---------------------------------------------------
  ComPWA::Tools::Plotting::DalitzPlot pl(kin, "DalitzFit", 100);
  pl.fillData(sample);
  pl.fillPhaseSpaceData(
      phspSample,
      std::make_shared<ComPWA::FunctionTree::FunctionTreeIntensity>(intens),
      "jpsiGammaPiPi", "", kBlue - 4);
  pl.plot();
  LOG(INFO) << "Done";

  return 0;
}
