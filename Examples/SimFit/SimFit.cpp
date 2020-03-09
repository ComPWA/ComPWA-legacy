// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/FunctionTree/FunctionTreeEstimator.hpp"
#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootGenerator.hpp"
#include "Physics/BuilderXML.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Tools/Plotting/RootPlotData.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Estimator/MinLogLH/SumMinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

std::string MyParticleList = R"####(
<ParticleList>
  <Particle Name="f0(980)">
    <Pid>9000111</Pid>
    <Parameter Type="Mass" Name="Mass_f0(980)">
      <Value>0.994</Value>
      <Error>0.001</Error>
      <Fix>false</Fix>
    </Parameter>
    <QuantumNumber Class="Spin" Type="Spin" Value="0"/>
    <QuantumNumber Class="Int" Type="Charge" Value="0"/>
    <QuantumNumber Class="Int" Type="Parity" Value="1"/>
    <DecayInfo Type="flatte">
      <FormFactor Type="0" />
      <Parameter Type="Coupling" Name="gPiPi_f0(980)">
        <Value>2.66</Value>
        <Error>0.001</Error>
        <Fix>true</Fix>
        <ParticleA>pi+</ParticleA>
        <ParticleB>pi-</ParticleB>
      </Parameter>
      <Parameter Type="Coupling" Name="gKK_f0(980)">
        <Value>3.121343843602647</Value>
        <Error>0.001</Error>
        <Fix>true</Fix>
        <ParticleA>K+</ParticleA>
        <ParticleB>K-</ParticleB>
      </Parameter>
      <Parameter Type="Coupling" Name="gKK_f0(980)">
        <Value>3.121343843602647</Value>
        <Error>0.001</Error>
        <Fix>true</Fix>
        <ParticleA>K_S0</ParticleA>
        <ParticleB>K_S0</ParticleB>
      </Parameter>
      <Parameter Type="MesonRadius" Name="Radius_f0(980)">
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name="Zc(3900)">
    <Pid>999999999</Pid>
    <Parameter Type="Mass" Name="Mass_Zc(3900)">
      <Value>3.900</Value>
      <Error>0.001</Error>
      <Fix>false</Fix>
    </Parameter>
    <QuantumNumber Class="Spin" Type="Spin" Value="0"/>
    <QuantumNumber Class="Int" Type="Charge" Value="0"/>
    <QuantumNumber Class="Int" Type="Parity" Value="1"/>
    <DecayInfo Type="flatte">
      <FormFactor Type="0" />
      <Parameter Type="Coupling" Name="gPiJPsi_Zc(3900)">
        <Value>0.075</Value>
        <Error>0.001</Error>
        <Fix>true</Fix>
        <ParticleA>pi+</ParticleA>
        <ParticleB>J/psi</ParticleB>
      </Parameter>
      <Parameter Type="Coupling" Name="gDD_Zc(3900)">
        <Value>2.03</Value>
        <Error>0.001</Error>
        <Fix>true</Fix>
        <ParticleA>D+</ParticleA>
        <ParticleB>D-</ParticleB>
      </Parameter>
      <Parameter Type="MesonRadius" Name="Radius_Zc(3900)">
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
</ParticleList>
)####";

std::string ModelSqrtS4230 = R"####(
<Intensity Class="IncoherentIntensity">
  <Intensity Class="CoherentIntensity">
    <Amplitude Class="CoefficientAmplitude">
      <Parameter Type="Magnitude" Name="Magnitude_f0(980)0">
        <Value>1.</Value>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Type="Phase" Name="Phase_f0(980)0">
        <Value>0.</Value>
        <Fix>true</Fix>
      </Parameter>
      <Amplitude Class="SequentialAmplitude">
        <Amplitude Class="HelicityDecay">
          <DecayParticle Name="f0(980)" Helicity="0"/>
          <RecoilSystem FinalState="2" />
          <DecayProducts>
            <Particle Name="pi+" FinalState="0"  Helicity="0"/>
            <Particle Name="pi-" FinalState="1"  Helicity="0"/>
          </DecayProducts>
        </Amplitude>
      </Amplitude>
    </Amplitude>
    <Amplitude Class="SequentialAmplitude">
      <Amplitude Class="HelicityDecay">
        <DecayParticle Name="Zc(3900)" Helicity="0"/>
        <RecoilSystem FinalState="0" />
        <DecayProducts>
          <Particle Name="pi-" FinalState="1"  Helicity="0"/>
          <Particle Name="J/psi" FinalState="2"  Helicity="0"/>
        </DecayProducts>
      </Amplitude>
    </Amplitude>
    <Amplitude Class="SequentialAmplitude">
      <Amplitude Class="HelicityDecay">
        <DecayParticle Name="Zc(3900)" Helicity="0"/>
        <RecoilSystem FinalState="1" />
        <DecayProducts>
          <Particle Name="pi+" FinalState="0"  Helicity="0"/>
          <Particle Name="J/psi" FinalState="2"  Helicity="0"/>
        </DecayProducts>
      </Amplitude>
    </Amplitude>
  </Intensity>
</Intensity>
)####";

std::string ModelSqrtS4260 = R"####(
<Intensity Class="IncoherentIntensity">
  <Intensity Class="CoherentIntensity">
    <Amplitude Class="CoefficientAmplitude">
      <Parameter Type="Magnitude" Name="Magnitude_f0(980)0">
        <Value>1.</Value>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Type="Phase" Name="Phase_f0(980)0">
        <Value>0.</Value>
        <Fix>true</Fix>
      </Parameter>
      <Amplitude Class="SequentialAmplitude">
        <Amplitude Class="HelicityDecay">
          <DecayParticle Name="f0(980)" Helicity="0"/>
          <RecoilSystem FinalState="2" />
          <DecayProducts>
            <Particle Name="pi+" FinalState="0"  Helicity="0"/>
            <Particle Name="pi-" FinalState="1"  Helicity="0"/>
          </DecayProducts>
        </Amplitude>
      </Amplitude>
    </Amplitude>
  </Intensity>
</Intensity>
)####";

using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

struct EnergyParameters {
  size_t NumberOfEvents;
  ComPWA::Physics::HelicityFormalism::HelicityKinematics Kinematics;
  ComPWA::FunctionTree::FunctionTreeIntensity Intensity;
  ComPWA::EventCollection DataSample;
  ComPWA::EventCollection PhspSample;
  ComPWA::Data::DataSet Points;
};

EnergyParameters createEnergyPar(size_t NumberOfEvents, double SqrtS, int Seed,
                                 ComPWA::ParticleList ParticleList,
                                 std::vector<ComPWA::pid> InitialState,
                                 std::vector<ComPWA::pid> FinalState,
                                 const boost::property_tree::ptree &ModelTree) {
  auto Kinematics = ComPWA::Physics::HelicityFormalism::HelicityKinematics(
      ParticleList, InitialState, FinalState,
      ComPWA::FourMomentum(0.0, 0.0, 0.0, SqrtS));

  ComPWA::Data::Root::RootGenerator Generator(
      Kinematics.getParticleStateTransitionKinematicsInfo());

  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(Seed);

  auto PhspSample =
      ComPWA::Data::generatePhsp(100000, Generator, RandomGenerator);

  // new builder because we need to use a different phsp sample for
  // normalization
  auto Builder = ComPWA::Physics::IntensityBuilderXML(
      ParticleList, Kinematics, ModelTree.get_child("Intensity"), PhspSample);
  // Construct intensity class from model string
  auto Intensity = Builder.createIntensity();

  auto DataSample = ComPWA::Data::generate(
      NumberOfEvents, Kinematics, RandomGenerator, Intensity, PhspSample);
  auto DataSet = Kinematics.convert(DataSample);
  return EnergyParameters{NumberOfEvents,       std::move(Kinematics),
                          std::move(Intensity), DataSample,
                          PhspSample,           DataSet};
}

///
/// Simultaneous fit of multiple energy points of the reaction
/// \f$ e^+e^- \to pi^+ pi^- J/psi \f$.
///
int main(int argc, char **argv) {
  using namespace ComPWA;
  using namespace boost::property_tree;

  // initialize logging
  Logging Log("log.txt", "debug");

  std::stringstream ModelStream;

  // List with all particle information needed
  auto ParticleList = readParticles("particle_list.xml");
  ModelStream << MyParticleList;
  insertParticles(ParticleList, ModelStream);

  std::vector<pid> InitialState({11, -11});
  std::vector<pid> FinalState({211, -211, 443});
  FourMomentum CmsP4;

  //---------------------------------------------------
  // sqrtS = 4230
  //---------------------------------------------------
  ModelStream.clear();
  ModelStream << ModelSqrtS4230;
  ptree ModelTr;
  xml_parser::read_xml(ModelStream, ModelTr);

  EnergyParameters SqrtS4230 = createEnergyPar(
      1000, 4.230, 123, ParticleList, InitialState, FinalState, ModelTr);

  //---------------------------------------------------
  // sqrtS = 4260
  //---------------------------------------------------
  ModelStream.clear();
  ModelStream << ModelSqrtS4260;
  ptree ModelTree4260;
  xml_parser::read_xml(ModelStream, ModelTree4260);

  EnergyParameters SqrtS4260 = createEnergyPar(
      1000, 4.260, 456, ParticleList, InitialState, FinalState, ModelTree4260);

  //---------------------------------------------------
  // Construct combined likelihood
  //---------------------------------------------------
  std::vector<
      std::pair<ComPWA::FunctionTree::FunctionTreeEstimator, FitParameterList>>
      Estimators;
  Estimators.push_back(ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
      SqrtS4230.Intensity, SqrtS4230.Points));
  Estimators.push_back(ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
      SqrtS4260.Intensity, SqrtS4260.Points));

  auto LogLHSumEstimator =
      ComPWA::Estimator::createSumMinLogLHFunctionTreeEstimator(
          std::move(Estimators));

  //---------------------------------------------------
  // Run fit
  //---------------------------------------------------

  using ComPWA::Optimizer::Minuit2::MinuitResult;

  //  std::get<0>(LogLHSumEstimator)->print(25);
  //  LOG(INFO) << "Fit parameter list: " <<
  //  std::get<1>(LogLHSumEstimator).to_str();
  auto Optimizer = Optimizer::Minuit2::MinuitIF();
  auto FitResult = Optimizer.optimize(std::get<0>(LogLHSumEstimator),
                                      std::get<1>(LogLHSumEstimator));

  LOG(INFO) << FitResult;

  //---------------------------------------------------
  // 5.1) Save the fit result
  //---------------------------------------------------
  std::ofstream OutputStream("fitResult.xml");
  boost::archive::xml_oarchive Archive(OutputStream);
  Archive << BOOST_SERIALIZATION_NVP(FitResult);

  //---------------------------------------------------
  // 6) Plot data sample and intensity
  //---------------------------------------------------
  ComPWA::Tools::Plotting::RootPlotData PlotData4230(
      SqrtS4230.Kinematics.getParticleStateTransitionKinematicsInfo(),
      "plot.root", "RECREATE");
  PlotData4230.createDirectory("sqrtS4230");
  PlotData4230.writeData(SqrtS4230.Kinematics.convert(SqrtS4230.DataSample));
  PlotData4230.writeIntensityWeightedPhspSample(
      SqrtS4230.Kinematics.convert(SqrtS4230.PhspSample), SqrtS4230.Intensity);

  ComPWA::Tools::Plotting::RootPlotData PlotData4260(
      SqrtS4260.Kinematics.getParticleStateTransitionKinematicsInfo(),
      "plot.root", "UPDATE");
  PlotData4260.createDirectory("sqrtS4230");
  PlotData4260.writeData(SqrtS4260.Kinematics.convert(SqrtS4260.DataSample));
  PlotData4260.writeIntensityWeightedPhspSample(
      SqrtS4260.Kinematics.convert(SqrtS4260.PhspSample), SqrtS4260.Intensity);

  LOG(INFO) << "Done";

  return 0;
}
