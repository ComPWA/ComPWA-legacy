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
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/BuilderXML.hpp"
#include "Physics/ParticleList.hpp"
#include "Tools/Plotting/RootPlotData.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Estimator/MinLogLH/SumMinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

std::string partList = R"####(
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

std::string modelSqrtS4230 = R"####(
<Intensity Class="IncoherentIntensity" Name="sqrtS4230_inc">
  <Intensity Class="CoherentIntensity" Name="sqrtS4230">
    <Amplitude Class="CoefficientAmplitude" Name="f0(980)">
      <Parameter Type="Magnitude" Name="Magnitude_f0(980)0">
        <Value>1.</Value>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Type="Phase" Name="Phase_f0(980)0">
        <Value>0.</Value>
        <Fix>true</Fix>
      </Parameter>
      <Amplitude Class="SequentialAmplitude" Name="f0(980)">
        <Amplitude Class="HelicityDecay" Name="f0(980)ToKK">
          <DecayParticle Name="f0(980)" Helicity="0"/>
          <RecoilSystem FinalState="2" />
          <DecayProducts>
            <Particle Name="pi+" FinalState="0"  Helicity="0"/>
            <Particle Name="pi-" FinalState="1"  Helicity="0"/>
          </DecayProducts>
        </Amplitude>
      </Amplitude>
    </Amplitude>
    <Amplitude Class="SequentialAmplitude" Name="Zc(3900)_JpsiPiMinus">
      <Amplitude Class="HelicityDecay" Name="Zc(3900)_JpsiPiMinusRes">
        <DecayParticle Name="Zc(3900)" Helicity="0"/>
        <RecoilSystem FinalState="0" />
        <DecayProducts>
          <Particle Name="pi-" FinalState="1"  Helicity="0"/>
          <Particle Name="J/psi" FinalState="2"  Helicity="0"/>
        </DecayProducts>
      </Amplitude>
    </Amplitude>
    <Amplitude Class="SequentialAmplitude" Name="Zc(3900)_JpsiPiPlus">
      <Amplitude Class="HelicityDecay" Name="Zc(3900)_JpsiPiPlusRes">
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

std::string modelSqrtS4260 = R"####(
<Intensity Class="IncoherentIntensity" Name="sqrtS4260_inc">
  <Intensity Class="CoherentIntensity" Name="sqrtS4260">
    <Amplitude Class="CoefficientAmplitude" Name="f0(980)">
      <Parameter Type="Magnitude" Name="Magnitude_f0(980)0">
        <Value>1.</Value>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Type="Phase" Name="Phase_f0(980)0">
        <Value>0.</Value>
        <Fix>true</Fix>
      </Parameter>
      <Amplitude Class="SequentialAmplitude" Name="f0(980)">
        <Amplitude Class="HelicityDecay"  Name="f0(980)ToKK">
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

struct energyPar {
  size_t _nEvents;
  ComPWA::Physics::HelicityFormalism::HelicityKinematics Kinematics;
  ComPWA::FunctionTree::FunctionTreeIntensity Intens;
  std::vector<ComPWA::Event> DataSample;
  std::vector<ComPWA::Event> McSample;
  ComPWA::Data::DataSet Points;
};

energyPar createEnergyPar(size_t NEvents, double SqrtS, int seed,
                          std::shared_ptr<ComPWA::PartList> partL,
                          std::vector<ComPWA::pid> initialState,
                          std::vector<ComPWA::pid> finalState,
                          const boost::property_tree::ptree &ModelPT) {
  auto kin = ComPWA::Physics::HelicityFormalism::HelicityKinematics(
      partL, initialState, finalState,
      ComPWA::FourMomentum(0.0, 0.0, 0.0, SqrtS));

  ComPWA::Data::Root::RootGenerator gen(
      kin.getParticleStateTransitionKinematicsInfo());

  ComPWA::Data::Root::RootUniformRealGenerator RandomGenerator(seed);

  auto mcSample = ComPWA::Data::generatePhsp(100000, gen, RandomGenerator);

  // new builder because we need to use a different phsp sample for
  // normalization
  auto Builder = ComPWA::Physics::IntensityBuilderXML(
      partL, kin, ModelPT.get_child("Intensity"), mcSample);
  // Construct intensity class from model string
  auto amp = Builder.createIntensity();

  auto data =
      ComPWA::Data::generate(NEvents, kin, RandomGenerator, amp, mcSample);
  auto dataSet = ComPWA::Data::convertEventsToDataSet(data, kin);
  return energyPar{NEvents, std::move(kin), std::move(amp),
                   data,    mcSample,       dataSet};
}

///
/// Simulaneous fit of multiple energy points of the reaction
/// e+e- \to pi+ pi - J/psi.
///
int main(int argc, char **argv) {
  using namespace ComPWA;
  using namespace boost::property_tree;

  // initialize logging
  Logging log("log.txt", "debug");

  ptree tmpTr;
  std::stringstream modelStream;

  // List with all particle information needed
  auto partL = std::make_shared<PartList>();
  ReadParticles(partL, ComPWA::Physics::defaultParticleList);
  modelStream << partList;
  xml_parser::read_xml(modelStream, tmpTr);
  modelStream.clear();
  ReadParticles(partL, tmpTr);

  std::vector<pid> initialState({11, -11});
  std::vector<pid> finalState({211, -211, 443});
  FourMomentum cmsP4;

  //---------------------------------------------------
  // sqrtS = 4230
  //---------------------------------------------------
  modelStream.clear();
  modelStream << modelSqrtS4230;
  ptree ModelTr;
  xml_parser::read_xml(modelStream, ModelTr);

  energyPar sqrtS4230 = createEnergyPar(1000, 4.230, 123, partL, initialState,
                                        finalState, ModelTr);

  //---------------------------------------------------
  // sqrtS = 4260
  //---------------------------------------------------
  modelStream.clear();
  modelStream << modelSqrtS4260;
  ptree ModelTr4260;
  xml_parser::read_xml(modelStream, ModelTr4260);

  energyPar sqrtS4260 = createEnergyPar(1000, 4.260, 456, partL, initialState,
                                        finalState, ModelTr4260);

  //---------------------------------------------------
  // Construct combined likelihood
  //---------------------------------------------------
  std::vector<
      std::pair<ComPWA::FunctionTree::FunctionTreeEstimator, FitParameterList>>
      v;
  v.push_back(ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
      sqrtS4230.Intens, sqrtS4230.Points));
  v.push_back(ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
      sqrtS4260.Intens, sqrtS4260.Points));

  auto LogLHSumEstimator =
      ComPWA::Estimator::createSumMinLogLHFunctionTreeEstimator(std::move(v));

  //---------------------------------------------------
  // Run fit
  //---------------------------------------------------

  using ComPWA::Optimizer::Minuit2::MinuitResult;

  //  std::get<0>(LogLHSumEstimator)->print(25);
  //  LOG(INFO) << "Fit parameter list: " <<
  //  std::get<1>(LogLHSumEstimator).to_str();
  auto minuitif = Optimizer::Minuit2::MinuitIF();
  auto result = minuitif.optimize(std::get<0>(LogLHSumEstimator),
                                  std::get<1>(LogLHSumEstimator));

  LOG(INFO) << result;

  //---------------------------------------------------
  // 5.1) Save the fit result
  //---------------------------------------------------
  std::ofstream ofs("fitResult.xml");
  boost::archive::xml_oarchive oa(ofs);
  oa << BOOST_SERIALIZATION_NVP(result);

  //---------------------------------------------------
  // 6) Plot data sample and intensity
  //---------------------------------------------------
  ComPWA::Tools::Plotting::RootPlotData sqrtS4230_pl(
      sqrtS4230.Kinematics.getParticleStateTransitionKinematicsInfo(),
      "plot.root", "RECREATE");
  sqrtS4230_pl.createDirectory("sqrtS4230");
  sqrtS4230_pl.writeData(
      Data::convertEventsToDataSet(sqrtS4230.DataSample, sqrtS4230.Kinematics));
  sqrtS4230_pl.writeIntensityWeightedPhspSample(
      Data::convertEventsToDataSet(sqrtS4230.McSample, sqrtS4230.Kinematics),
      sqrtS4230.Intens);

  ComPWA::Tools::Plotting::RootPlotData sqrtS4260_pl(
      sqrtS4260.Kinematics.getParticleStateTransitionKinematicsInfo(),
      "plot.root", "UPDATE");
  sqrtS4260_pl.createDirectory("sqrtS4230");
  sqrtS4260_pl.writeData(
      Data::convertEventsToDataSet(sqrtS4260.DataSample, sqrtS4260.Kinematics));
  sqrtS4260_pl.writeIntensityWeightedPhspSample(
      Data::convertEventsToDataSet(sqrtS4260.McSample, sqrtS4260.Kinematics),
      sqrtS4260.Intens);

  LOG(INFO) << "Done";

  return 0;
}
