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
#include "Physics/BuilderXML.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
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
std::string AmplitudeModel = R"####(
<Intensity Class="NormalizedIntensity">
  <IntegrationStrategy Class="MCIntegrationStrategy"/>
  <Intensity Class="CoherentIntensity" Component="jpsiGammaPiPi">
    <Amplitude Class="CoefficientAmplitude" Component="f2(1270)">
	    <Parameter Class='Double' Type="Magnitude" Name="Magnitude_f2">
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
      <Amplitude Class="NormalizedAmplitude">
        <IntegrationStrategy Class="MCIntegrationStrategy"/>
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
	  <Amplitude Class="CoefficientAmplitude" Component="myAmp">
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
      <Amplitude Class="NormalizedAmplitude">
        <IntegrationStrategy Class="MCIntegrationStrategy"/>
		    <Amplitude Class="HelicityDecay" >
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
	<Parameter Type="Mass" Name="Mass_jpsi">
	  <Value>3.096900</Value>
	  <Fix>true</Fix>
	</Parameter>
	<QuantumNumber Class="Spin" Type="Spin" Value="1" />
	<QuantumNumber Class="Int" Type="Charge" Value="0" />
	<QuantumNumber Class="Int" Type="Parity" Value="-1" />
	<QuantumNumber Class="Int" Type="Cparity" Value="-1" />
	<QuantumNumber Class="Int" Type="Gparity" Value="-1" />
	<QuantumNumber Class="Spin" Type="IsoSpin" Value="0" Projection="0" />
	<QuantumNumber Class="Int" Type="BaryonNumber" Value="0" />
	<QuantumNumber Class="Int" Type="Charm" Value="0" />
	<QuantumNumber Class="Int" Type="Strangeness" Value="0" />
	<DecayInfo Type="relativisticBreitWigner">
	  <FormFactor Type="0" />
	  <Parameter Type="Width" Name="Width_jpsi">
		<Value>9.29E-05</Value>
		<Fix>true</Fix>
	  </Parameter>
	  <Parameter Type="MesonRadius" Name="Radius_jpsi">
		<Value>2.5</Value>
		<Fix>true</Fix>
		<Min>2.0</Min>
		<Max>3.0</Max>
	  </Parameter>
	</DecayInfo>
  </Particle>
  <Particle Name="pi0">
	<Pid>111</Pid>
	<Parameter Type="Mass" Name="Mass_neutralPion">
	  <Value>0.1349766</Value>
	  <Error>0.000006</Error>
	</Parameter>
	<QuantumNumber Class="Spin" Type="Spin" Value="0" />
	<QuantumNumber Class="Int" Type="Charge" Value="0" />
	<QuantumNumber Class="Int" Type="Parity" Value="-1" />
	<QuantumNumber Class="Int" Type="Cparity" Value="1" />
	<QuantumNumber Class="Int" Type="Gparity" Value="-1" />
	<QuantumNumber Class="Spin" Type="IsoSpin" Value="1" Projection="0" />
	<QuantumNumber Class="Int" Type="BaryonNumber" Value="0" />
	<QuantumNumber Class="Int" Type="Charm" Value="0" />
	<QuantumNumber Class="Int" Type="Strangeness" Value="0" />
  </Particle>
  <Particle Name="gamma">
	<Pid>22</Pid>
	<Parameter Type="Mass" Name="Mass_gamma">
	  <Value>0.0</Value>
	</Parameter>
	<QuantumNumber Class="Spin" Type="Spin" Value="1" />
	<QuantumNumber Class="Int" Type="Charge" Value="0" />
	<QuantumNumber Class="Int" Type="Parity" Value="-1" />
	<QuantumNumber Class="Int" Type="Cparity" Value="-1" />
	<QuantumNumber Class="Spin" Type="IsoSpin" Value="0" Projection="0" />
	<QuantumNumber Class="Int" Type="BaryonNumber" Value="0" />
	<QuantumNumber Class="Int" Type="Charm" Value="0" />
	<QuantumNumber Class="Int" Type="Strangeness" Value="0" />
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
  Logging Log("debug", "DalitzFit-log.txt");

  // List with all particle information needed
  std::stringstream ParticlesStream(MyParticleList);
  ParticleList Particles = readParticles(ParticlesStream);

  //---------------------------------------------------
  // 1) Create Kinematics object
  //---------------------------------------------------
  std::vector<pid> InitialState = {443};
  std::vector<pid> FinalState = {22, 111, 111};
  HelicityKinematics Kinematics(Particles, InitialState, FinalState);

  //---------------------------------------------------
  // 2) Generate a large phase space sample
  //---------------------------------------------------
  ComPWA::Data::Root::RootGenerator Generator(
      Kinematics.getParticleStateTransitionKinematicsInfo());

  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(173);

  auto PhspSample(Data::generatePhsp(100000, Generator, RandomGenerator));

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
      Particles, Kinematics, ModelTree.get_child("Intensity"), PhspSample);
  auto Intensity = Builder.createIntensity();

  //---------------------------------------------------
  // 4) Generate a data sample given intensity and kinematics
  //---------------------------------------------------
  auto DataSample = ComPWA::Data::generate(1000, Kinematics, RandomGenerator,
                                           Intensity, PhspSample);

  auto SampleDataSet = Kinematics.convert(DataSample);
  //---------------------------------------------------
  // 5) Fit the model to the data and print the result
  //---------------------------------------------------
  auto Estimator = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
      Intensity, SampleDataSet);

  LOG(DEBUG) << Estimator.first.print(25);

  auto Optimizer = Optimizer::Minuit2::MinuitIF();

  // STARTING MINIMIZATION
  auto FitResult =
      Optimizer.optimize(std::get<0>(Estimator), std::get<1>(Estimator));

  LOG(INFO) << FitResult;

  // calculate fit fractions and errors
  auto Components = Builder.createIntensityComponents(
      {{"f2(1270)"}, {"myAmp"}, {"jpsiGammaPiPi"}});

  auto MyFractions = {std::make_pair(Components[0], Components[2]),
                      std::make_pair(Components[1], Components[2])};

  ComPWA::Tools::FitFractions FitFraction;
  auto FitFractionList =
      FitFraction.calculateFitFractionsWithCovarianceErrorPropagation(
          MyFractions, Kinematics.convert(PhspSample), FitResult);

  LOG(INFO) << FitFractionList;

  //---------------------------------------------------
  // 5.1) Save the fit result
  //---------------------------------------------------
  std::ofstream FileStream("DalitzFit-fitResult.xml");
  boost::archive::xml_oarchive Archive(FileStream);
  Archive << BOOST_SERIALIZATION_NVP(FitResult);

  //---------------------------------------------------
  // 6) Plot data sample and intensity
  //---------------------------------------------------
  ComPWA::Tools::Plotting::DalitzPlot Plot(Kinematics, "DalitzFit", 100);
  Plot.fill(DataSample, true, "data", "Data sample", kBlack);
  Plot.fill(PhspSample, false, "phsp", "Phsp sample", kGreen);
  Plot.fill(PhspSample, Intensity, false, "fit", "jpsiGammaPiPi model", kBlue);
  Plot.plot();
  LOG(INFO) << "Done";

  return 0;
}
