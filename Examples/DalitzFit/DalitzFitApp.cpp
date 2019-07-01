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

#include "Core/FunctionTreeIntensityWrapper.hpp"
#include "Core/Intensity.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Physics/CoherentIntensity.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Physics/NormalizationIntensityDecorator.hpp"
#include "Physics/ParticleList.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Generate.hpp"
#include "Tools/ParameterTools.hpp"
#include "Tools/Plotting/DalitzPlot.hpp"
#include "Tools/RootGenerator.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

using namespace ComPWA;
using ComPWA::Optimizer::Minuit2::MinuitResult;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

// Enable serialization of MinuitResult. For some reason has to be outside
// any namespaces.
BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)

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
  Logging log("DalitzFit-log.txt", "debug");

  // List with all particle information needed
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, ComPWA::Physics::defaultParticleList);
  ReadParticles(partL, myParticles);

  //---------------------------------------------------
  // 1) Create Kinematics object
  //---------------------------------------------------
  std::vector<pid> initialState = {443};
  std::vector<pid> finalState = {22, 111, 111};
  auto kin =
      std::make_shared<HelicityKinematics>(partL, initialState, finalState);

  //---------------------------------------------------
  // 2) Generate a large phase space sample
  //---------------------------------------------------
  auto gen = std::make_shared<ComPWA::Tools::RootGenerator>(
      kin->getParticleStateTransitionKinematicsInfo(), 173);
  std::shared_ptr<ComPWA::Data::DataSet> phspSample(
      ComPWA::Tools::generatePhsp(1000000, gen));

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
      Builder.createOldIntensity(partL, kin, modelTree.get_child("Intensity"));

  // Pass phsp sample to intensity for normalization.
  // Convert to dataPoints first.
  // auto phspPoints =
  //    std::make_shared<std::vector<DataPoint>>(phspSample->dataPoints(kin));

  auto newIntens =
      std::make_shared<ComPWA::FunctionTreeIntensityWrapper>(intens, kin);

  //---------------------------------------------------
  // 4) Generate a data sample given intensity and kinematics
  //---------------------------------------------------
  std::shared_ptr<ComPWA::Data::DataSet> sample =
      ComPWA::Tools::generate(1000, kin, gen, newIntens, phspSample);

  sample->convertEventsToParameterList(kin);
  phspSample->convertEventsToParameterList(kin);

  //---------------------------------------------------
  // 5) Fit the model to the data and print the result
  //---------------------------------------------------
  ParameterList fitPar;
  intens->addUniqueParametersTo(fitPar);
  // Set start error of 0.05 for parameters, run Minos?
  setErrorOnParameterList(fitPar, 0.05, false);

  auto esti =
      ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(intens, sample);

  // LOG(DEBUG) << esti->print(25);

  auto minuitif = new Optimizer::Minuit2::MinuitIF(esti, fitPar);
  minuitif->setUseHesse(true);

  // STARTING MINIMIZATION
  auto result = std::dynamic_pointer_cast<MinuitResult>(minuitif->exec(fitPar));

  // First extract the normalized intensity, and then the coherent intensity.
  try {
    // before we calculate the fit fractions we have to update the fit
    // parameters of the intensity (since we used a function tree for fitting)
    intens->updateParametersFrom(result->finalParameters());
    auto NormIntens = std::dynamic_pointer_cast<
        ComPWA::Physics::NormalizationIntensityDecorator>(intens);
    auto CohIntens =
        std::dynamic_pointer_cast<ComPWA::Physics::CoherentIntensity>(
            NormIntens->getUnnormalizedIntensity());

    std::vector<std::string> AmpComponents = {"myAmp", "f2(1270)"};
    ParameterList fitFracs =
        Tools::calculateFitFractionsWithCovarianceErrorPropagation(
            CohIntens, phspSample, result->covarianceMatrix(), AmpComponents);
    result->setFitFractions(fitFracs);
  } catch (const std::exception &e) {
    throw std::runtime_error(e.what());
  }
  result->print();

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
  pl.setData(sample);
  pl.setPhspData(phspSample);
  pl.setFitAmp(newIntens, "jpsiGammaPiPi", "", kBlue - 4);
  pl.plot();
  LOG(INFO) << "Done";

  return 0;
}
