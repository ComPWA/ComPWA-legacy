// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

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

#include "Core/FunctionTreeEstimatorWrapper.hpp"
#include "Core/FunctionTreeIntensityWrapper.hpp"
#include "Core/Intensity.hpp"
#include "Core/Logging.hpp"
#include "Core/ProgressBar.hpp"
#include "Core/Properties.hpp"
#include "Data/DataSet.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Physics/ParticleList.hpp"
#include "Tools/Generate.hpp"
#include "Tools/ParameterTools.hpp"

#include "Tools/Plotting/RootPlotData.hpp"
#include "Tools/RootGenerator.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Estimator/MinLogLH/SumMinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

using namespace boost::property_tree;

using namespace ComPWA;
using namespace ComPWA::Tools;
using namespace ComPWA::Physics::HelicityFormalism;
using ComPWA::Optimizer::Minuit2::MinuitResult;

// Enable serialization of MinuitResult. For some reason has to be outside
// any namespaces.
BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)

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

struct energyPar {
  int _nEvents;
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::HelicityKinematics> _kin;
  std::shared_ptr<ComPWA::Generator> _gen;
  std::shared_ptr<ComPWA::OldIntensity> _amp;
  std::shared_ptr<ComPWA::Data::DataSet> _data;
  std::shared_ptr<ComPWA::Data::DataSet> _mcSample;
  std::shared_ptr<ComPWA::Tools::Plotting::RootPlotData> _pl;
};

///
/// Simulaneous fit of multiple energy points of the reaction
/// e+e- \to pi+ pi - J/psi.
///
int main(int argc, char **argv) {

  // initialize logging
  Logging log("log.txt", "debug");

  ptree tmpTr;
  std::stringstream modelStream;

  // List with all particle information needed
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, ComPWA::Physics::defaultParticleList);
  modelStream << partList;
  xml_parser::read_xml(modelStream, tmpTr);
  modelStream.clear();
  ReadParticles(partL, tmpTr);
  ParameterList fitPar;

  std::vector<pid> initialState({11, -11});
  std::vector<pid> finalState({211, -211, 443});
  FourMomentum cmsP4;

  //---------------------------------------------------
  // sqrtS = 4230
  //---------------------------------------------------
  energyPar sqrtS4230;
  cmsP4 = FourMomentum(0, 0, 0, 4.230);
  sqrtS4230._nEvents = 1000;

  sqrtS4230._kin = std::make_shared<HelicityKinematics>(partL, initialState,
                                                        finalState, cmsP4);
  sqrtS4230._gen = std::make_shared<ComPWA::Tools::RootGenerator>(
      sqrtS4230._kin->getParticleStateTransitionKinematicsInfo(), 123);

  modelStream << modelSqrtS4230;
  xml_parser::read_xml(modelStream, tmpTr);
  modelStream.clear();

  //  sqrtS4230._data =
  //      std::make_shared<ComPWA::DataReader::RootReader>("data4230.root",
  //      "data");
  sqrtS4230._mcSample = ComPWA::Tools::generatePhsp(100000, sqrtS4230._gen);

  // Construct intensity class from model string
  ComPWA::Physics::IntensityBuilderXML Builder;
  sqrtS4230._amp = Builder.createOldIntensity(partL, sqrtS4230._kin,
                                              tmpTr.get_child("Intensity"));
  sqrtS4230._amp->addUniqueParametersTo(fitPar);

  auto newAmp = std::make_shared<ComPWA::FunctionTreeIntensityWrapper>(
      sqrtS4230._amp, sqrtS4230._kin);

  sqrtS4230._data =
      ComPWA::Tools::generate(sqrtS4230._nEvents, sqrtS4230._kin,
                              sqrtS4230._gen, newAmp, sqrtS4230._mcSample);

  sqrtS4230._mcSample->convertEventsToParameterList(sqrtS4230._kin);
  auto estimator1 = ComPWA::Estimator::createMinLogLHEstimatorFunctionTree(
      sqrtS4230._amp, sqrtS4230._data, sqrtS4230._mcSample);
  // estimator1->head()->print();

  //---------------------------------------------------
  // sqrtS = 4260
  //---------------------------------------------------
  energyPar sqrtS4260;
  cmsP4 = FourMomentum(0, 0, 0, 4.260);
  sqrtS4260._nEvents = 1000;

  sqrtS4260._kin = std::make_shared<HelicityKinematics>(partL, initialState,
                                                        finalState, cmsP4);
  sqrtS4260._gen = std::make_shared<ComPWA::Tools::RootGenerator>(
      sqrtS4260._kin->getParticleStateTransitionKinematicsInfo(), 456);

  modelStream << modelSqrtS4260;
  xml_parser::read_xml(modelStream, tmpTr);

  //  sqrtS4260._data =
  //      std::make_shared<ComPWA::DataReader::RootReader>("data4230.root",
  //      "data");
  sqrtS4260._mcSample = ComPWA::Tools::generatePhsp(100000, sqrtS4260._gen);

  // Construct intensity class from model string
  sqrtS4260._amp = Builder.createOldIntensity(partL, sqrtS4260._kin,
                                              tmpTr.get_child("Intensity"));
  sqrtS4260._amp->addUniqueParametersTo(fitPar);

  auto newAmp2 = std::make_shared<ComPWA::FunctionTreeIntensityWrapper>(
      sqrtS4260._amp, sqrtS4260._kin);

  sqrtS4260._data =
      ComPWA::Tools::generate(sqrtS4260._nEvents, sqrtS4260._kin,
                              sqrtS4260._gen, newAmp2, sqrtS4260._mcSample);

  sqrtS4260._mcSample->convertEventsToParameterList(sqrtS4260._kin);
  auto estimator2 = ComPWA::Estimator::createMinLogLHEstimatorFunctionTree(
      sqrtS4260._amp, sqrtS4260._data, sqrtS4260._mcSample);
  // estimator2->head()->print();

  //---------------------------------------------------
  // sqrtS = 4340
  //---------------------------------------------------

  auto LogLHSumEstimator =
      std::make_shared<ComPWA::FunctionTreeEstimatorWrapper>(
          ComPWA::Estimator::createSumMinLogLHEstimatorFunctionTree(
              {estimator1, estimator2}),
          fitPar);

  //---------------------------------------------------
  // Run fit
  //---------------------------------------------------

  // LogLHSumEstimator->print(25);
  LOG(INFO) << "Fit parameter list: " << fitPar.to_str();
  auto minuitif = new Optimizer::Minuit2::MinuitIF(LogLHSumEstimator, fitPar);
  minuitif->setUseHesse(true);
  auto result = minuitif->exec(fitPar);

  result->print();

  //---------------------------------------------------
  // 5.1) Save the fit result
  //---------------------------------------------------
  std::ofstream ofs("fitResult.xml");
  boost::archive::xml_oarchive oa(ofs);
  oa << BOOST_SERIALIZATION_NVP(result);

  UpdateParticleList(partL, fitPar);
  ptree ptout;
  ptout.add_child("ParticleList", SaveParticles(partL));
  xml_parser::write_xml("fitParticles.xml", ptout, std::locale());

  //---------------------------------------------------
  // 6) Plot data sample and intensity
  //---------------------------------------------------
  /*sqrtS4230._pl = std::make_shared<ComPWA::Tools::Plotting::RootPlotData>(
      sqrtS4230._kin, sqrtS4230._amp);
  sqrtS4230._pl->setData(sqrtS4230._data);
  sqrtS4230._pl->setPhspMC(sqrtS4230._mcSample);
  sqrtS4230._pl->addComponent("f0(980)", "f0_980");
  sqrtS4230._pl->addComponent("Zc(3900)_JpsiPiPlus + Zc(3900)_JpsiPiMinus",
                              "Zc3900");
  sqrtS4230._pl->write("sqrtS4230", "plot.root", "RECREATE");

  sqrtS4260._pl = std::make_shared<ComPWA::Tools::Plotting::RootPlotData>(
      sqrtS4260._kin, sqrtS4260._amp);
  sqrtS4260._pl->setData(sqrtS4260._data);
  sqrtS4260._pl->setPhspMC(sqrtS4260._mcSample);
  sqrtS4260._pl->addComponent("f0(980)", "f0_980");
  sqrtS4260._pl->write("sqrtS4260", "plot.root", "UPDATE");
*/
  LOG(INFO) << "Done";

  return 0;
}
