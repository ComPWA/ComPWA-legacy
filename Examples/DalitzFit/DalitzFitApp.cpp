

// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Simple Dalitz plot analysis with ComPWA
///

#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "Physics/ParticleList.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/RunManager.hpp"
#include "Tools/DalitzPlot.hpp"
#include "Tools/ParameterTools.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

using namespace ComPWA;
using namespace ComPWA::DataReader;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;
using ComPWA::Physics::HelicityFormalism::IncoherentIntensity;
using ComPWA::Optimizer::Minuit2::MinuitResult;

// Enable serialization of MinuitResult. For some reason has to be outside
// any namespaces.
BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)

// We define an intensity model using a raw string literal. Currently, this is
// just a toy model without any physical meaning.
// (comments within the string are ignored!). This is convenient since we
// do not have to configure the build system to copy input files somewhere.
// In practise you may want to use a normal XML input file instead.
std::string amplitudeModel = R"####(
<IncoherentIntensity Name="jpsiGammaPiPi_inc">
	<CoherentIntensity Name="jpsiGammaPiPi">
  	<Amplitude Name="f2(1270)">
			<Parameter Type="Magnitude"	Name="Magnitude_f2">
				<Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>2.0</Max>
				<Fix>false</Fix>
			</Parameter>
			<Parameter Type="Phase" Name="Phase_f2">
				<Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
				<Fix>false</Fix>
			</Parameter>
			<Resonance Name="f2ToPiPi">
				<DecayParticle Name="f2(1270)" Helicity="0"/>
				<SubSystem>
					<RecoilSystem FinalState="0" />
					<DecayProducts>
						<Particle Name="pi0" FinalState="1"  Helicity="0"/>
						<Particle Name="pi0" FinalState="2"  Helicity="0"/>
					</DecayProducts>
				</SubSystem>
			</Resonance>
		</Amplitude>
		<Amplitude Name="myAmp">
			<Parameter Type="Magnitude"	Name="Magnitude_my">
				<Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>2.0</Max>
				<Fix>true</Fix>
			</Parameter>
			<Parameter Type="Phase" Name="Phase_my`">
				<Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
				<Fix>true</Fix>
			</Parameter>
			<Resonance Name="MyResToPiPi">
				<DecayParticle Name="myRes" Helicity="0"/>
				<SubSystem>
					<RecoilSystem FinalState="0" />
					<DecayProducts>
						<Particle Name="pi0" FinalState="1"  Helicity="0"/>
						<Particle Name="pi0" FinalState="2"  Helicity="0"/>
					</DecayProducts>
				</SubSystem>
			</Resonance>
		</Amplitude>
	</CoherentIntensity>
</IncoherentIntensity>
)####";

std::string myParticles = R"####(
<ParticleList>
	<Particle Name="f2(1270)">
		<Pid>225</Pid>
		<Parameter Type="Mass" Name="Mass_f2(1270)">
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
			<Parameter Type="Width" Name="Width_f2(1270)">
				<Value>0.1867</Value>
			</Parameter>
			<Parameter Type="MesonRadius" Name="Radius_rho">
				<Value>2.5</Value>
				<Fix>true</Fix>
			</Parameter>
		</DecayInfo>
	</Particle>
	<Particle Name="myRes">
		<Pid>999999</Pid>
		<Parameter Type="Mass" Name="Mass_myRes">
			<Value>2.0</Value>
			<Error>8.0E-04</Error>
		</Parameter>
		<QuantumNumber Class="Spin" Type="Spin" Value="1"/>
		<QuantumNumber Class="Int" Type="Charge" Value="0"/>
		<QuantumNumber Class="Int" Type="Parity" Value="+1"/>
		<QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
		<DecayInfo Type="relativisticBreitWigner">
			<FormFactor Type="0" />
			<Parameter Type="Width" Name="Width_myRes">
				<Value>1.0</Value>
        <Min>0.1</Min>
        <Max>1.0</Max>
        <Fix>false</Fix>
			</Parameter>
			<Parameter Type="MesonRadius" Name="Radius_myRes">
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
  Logging log("DalitzFit-log.txt", boost::log::trivial::debug);

  // RunManager is supposed to manage all objects of the analysis. It generates
  // data and starts the fitting procedure.
  ComPWA::Tools::RunManager run;

  // List with all particle information needed
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, defaultParticleList);
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
  auto gen = std::make_shared<ComPWA::Tools::RootGenerator>(partL, kin);
  run.SetGenerator(gen);
  std::shared_ptr<Data> phspSample(new Data());
  run.SetPhspSample(phspSample);
  run.GeneratePhsp(100000);

  //---------------------------------------------------
  // 3) Create intensity from pre-defined model
  //---------------------------------------------------
  // Read in model property_tree
  std::stringstream modelStream;
  modelStream << amplitudeModel;
  boost::property_tree::ptree modelTree;
  boost::property_tree::xml_parser::read_xml(modelStream, modelTree);

  // Construct intensity class from model string
  auto intens = IncoherentIntensity::Factory(
      partL, kin, modelTree.get_child("IncoherentIntensity"));

  // Pass phsp sample to intensity for normalization.
  // Convert to dataPoints first.
  auto phspPoints =
      std::make_shared<std::vector<dataPoint>>(phspSample->GetDataPoints(kin));
  intens->SetPhspSample(phspPoints, phspPoints);
  run.SetAmplitude(intens);

  //---------------------------------------------------
  // 4) Generate a data sample given intensity and kinematics
  //---------------------------------------------------
  std::shared_ptr<Data> sample(new Data());
  run.SetData(sample);
  run.Generate(kin, 1000);

  //---------------------------------------------------
  // 5) Fit the model to the data and print the result
  //---------------------------------------------------
  ParameterList fitPar;
  intens->GetParameters(fitPar);
  // Set start error of 0.05 for parameters, run Minos?
  setErrorOnParameterList(fitPar, 0.05, false);

  auto esti = std::make_shared<Estimator::MinLogLH>(
      kin, intens, sample, phspSample, phspSample, 0, 0);

  esti->UseFunctionTree(true);
  LOG(debug) << esti->GetTree()->Head()->Print(25);

  auto minuitif = new Optimizer::Minuit2::MinuitIF(esti, fitPar);
  minuitif->SetHesse(true);
  run.SetOptimizer(std::shared_ptr<Optimizer::Optimizer>(minuitif));

  // STARTING MINIMIZATION
  auto result = std::dynamic_pointer_cast<MinuitResult>(run.Fit(fitPar));

  // Calculate fit fractions
  std::vector<std::pair<std::string, std::string>> fitComponents;
  fitComponents.push_back(
      std::pair<std::string, std::string>("myAmp", "jpsiGammaPiPi"));
  fitComponents.push_back(
      std::pair<std::string, std::string>("f2(1270)", "jpsiGammaPiPi"));

  ParameterList fitFracs =
      Tools::CalculateFitFractions(kin, intens, phspPoints, fitComponents);
  // A proper calculation of the fit fraction uncertainty requires
  // the uncertainty of the fit parameters propagated. We do this
  // using a numerical approach. Using the covariance matrix
  // 100 independend sets of fit parameters are generated and the fit fractions
  // are recalculated. In the end we take the RMS.
  Tools::CalcFractionError(fitPar, result->GetCovarianceMatrix(), fitFracs, kin,
                           intens, phspPoints, 100, fitComponents);

  result->SetFitFractions(fitFracs);
  result->Print();

  //---------------------------------------------------
  // 5.1) Save the fit result
  //---------------------------------------------------
  std::ofstream ofs("DalitzFit-fitResult.xml");
  boost::archive::xml_oarchive oa(ofs);
  oa << BOOST_SERIALIZATION_NVP(result);

  UpdateParticleList(partL, fitPar);
  boost::property_tree::ptree ptout;
  ptout.add_child("ParticleList", SaveParticles(partL));
  ptout.add_child("IncoherentIntensity", IncoherentIntensity::Save(intens));
  boost::property_tree::xml_parser::write_xml("DalitzFit-Model.xml", ptout,
                                              std::locale());
  //---------------------------------------------------
  // 6) Plot data sample and intensity
  //---------------------------------------------------
  ComPWA::Tools::DalitzPlot pl(kin, "DalitzFit", 100);
  pl.SetData(sample);
  pl.SetPhspData(phspSample);
  pl.SetFitAmp(intens, "", kBlue - 4);
  pl.Plot();
  LOG(info) << "Done";

  return 0;
}
