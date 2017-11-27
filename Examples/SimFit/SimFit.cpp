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
#include "Core/ProgressBar.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/ParticleList.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/Generate.hpp"
#include "Tools/DalitzPlot.hpp"
#include "Tools/RootPlot.hpp"
#include "Tools/ParameterTools.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Estimator/MinLogLH/SumMinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

using namespace boost::property_tree;

using namespace ComPWA;
using namespace ComPWA::DataReader;
using namespace ComPWA::Tools;
using namespace ComPWA::Physics::HelicityFormalism;
using ComPWA::Physics::HelicityFormalism::IncoherentIntensity;
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
<IncoherentIntensity Name="sqrtS4230_inc">
  <Parameter Type="Strength" Name="strength_sqrtS4230_inc">
  <Value>1.</Value>
  <Fix>true</Fix>
  </Parameter>
  <CoherentIntensity Name="sqrtS4230">
  <Parameter Type="Strength" Name="strength_d0tokkk">
    <Value>1.</Value>
    <Fix>true</Fix>
  </Parameter>
  <Amplitude Name="f0(980)">
    <Parameter Type="Magnitude" Name="Magnitude_f0(980)0">
    <Value>1.</Value>
    <Fix>true</Fix>
    </Parameter>
    <Parameter Type="Phase" Name="Phase_f0(980)0">
    <Value>0.</Value>
    <Fix>true</Fix>
    </Parameter>
    <Resonance Name="f0(980)ToKK">
    <DecayParticle Name="f0(980)" Helicity="0"/>
    <SubSystem>
      <RecoilSystem FinalState="2" />
      <DecayProducts>
      <Particle Name="pi+" FinalState="0"  Helicity="0"/>
      <Particle Name="pi-" FinalState="1"  Helicity="0"/>
      </DecayProducts>
    </SubSystem>
    </Resonance>
  </Amplitude>
  <Amplitude Name="Zc(3900)_JpsiPiMinus">
    <Parameter Type="Magnitude" Name="Magnitude_Zc(3900)_JpsiPiMinus">
    <Value>1.</Value>
    <Fix>true</Fix>
    </Parameter>
    <Parameter Type="Phase" Name="Phase_Zc(3900)_JpsiPiMinus">
    <Value>0.</Value>
    <Fix>true</Fix>
    </Parameter>
    <Resonance Name="Zc(3900)_JpsiPiMinusRes">
    <DecayParticle Name="Zc(3900)" Helicity="0"/>
    <SubSystem>
      <RecoilSystem FinalState="0" />
      <DecayProducts>
        <Particle Name="pi-" FinalState="1"  Helicity="0"/>
        <Particle Name="J/psi" FinalState="2"  Helicity="0"/>
      </DecayProducts>
    </SubSystem>
    </Resonance>
  </Amplitude>
  <Amplitude Name="Zc(3900)_JpsiPiPlus">
    <Parameter Type="Magnitude" Name="Magnitude_Zc(3900)_JpsiPiPlus">
    <Value>1.</Value>
    <Fix>true</Fix>
    </Parameter>
    <Parameter Type="Phase" Name="Phase_Zc(3900)_JpsiPiPlus">
    <Value>0.</Value>
    <Fix>true</Fix>
    </Parameter>
    <Resonance Name="Zc(3900)_JpsiPiPlusRes">
    <DecayParticle Name="Zc(3900)" Helicity="0"/>
    <SubSystem>
      <RecoilSystem FinalState="1" />
      <DecayProducts>
        <Particle Name="pi+" FinalState="0"  Helicity="0"/>
        <Particle Name="J/psi" FinalState="2"  Helicity="0"/>
      </DecayProducts>
    </SubSystem>
    </Resonance>
  </Amplitude>
  </CoherentIntensity>
</IncoherentIntensity>
)####";

std::string modelSqrtS4260 = R"####(
<IncoherentIntensity Name="sqrtS4260_inc">
  <Parameter Type="Strength" Name="strength_sqrtS4260_inc">
  <Value>1.</Value>
  <Fix>true</Fix>
  </Parameter>
  <CoherentIntensity Name="sqrtS4260">
  <Parameter Type="Strength" Name="strength_d0tokkk">
    <Value>1.</Value>
    <Fix>true</Fix>
  </Parameter>
  <Amplitude Name="f0(980)">
    <Parameter Type="Magnitude" Name="Magnitude_f0(980)0">
    <Value>1.</Value>
    <Fix>true</Fix>
    </Parameter>
    <Parameter Type="Phase" Name="Phase_f0(980)0">
    <Value>0.</Value>
    <Fix>true</Fix>
    </Parameter>
    <Resonance Name="f0(980)ToKK">
    <Parameter Type="Magnitude" Name="Magnitude_f0(980)ToPiPi">
      <Value>1.</Value>
      <Fix>true</Fix>
    </Parameter>
    <Parameter Type="Phase" Name="Phase_f0(980)ToPiPi">
      <Value>0.</Value>
      <Fix>true</Fix>
    </Parameter>
    <DecayParticle Name="f0(980)" Helicity="0"/>
    <SubSystem>
      <RecoilSystem FinalState="2" />
      <DecayProducts>
      <Particle Name="pi+" FinalState="0"  Helicity="0"/>
      <Particle Name="pi-" FinalState="1"  Helicity="0"/>
      </DecayProducts>
    </SubSystem>
    </Resonance>
  </Amplitude>
  </CoherentIntensity>
</IncoherentIntensity>
)####";

struct energyPar {
  int _nEvents;
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::HelicityKinematics> _kin;
  std::shared_ptr<ComPWA::Generator> _gen;
  std::shared_ptr<ComPWA::AmpIntensity> _amp;
  std::shared_ptr<ComPWA::DataReader::Data> _data;
  std::shared_ptr<ComPWA::DataReader::Data> _mcSample;
  std::shared_ptr<ComPWA::DataReader::Data> _mcSampleTrue;
  std::shared_ptr<std::vector<ComPWA::dataPoint>> _mcPoints;
  std::shared_ptr<ComPWA::Efficiency> _eff;
  std::shared_ptr<ComPWA::Estimator::MinLogLH> _minLH;
  std::shared_ptr<RootPlot> _pl;
};

///
/// Simulaneous fit of multiple energy points of the reaction
/// e+e- \to pi+ pi - J/psi.
///
int main(int argc, char **argv) {

  // initialize logging
  Logging log("log.txt", boost::log::trivial::debug);

  ptree tmpTr;
  std::stringstream modelStream;

  // List with all particle information needed
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, defaultParticleList);
  modelStream << partList;
  xml_parser::read_xml(modelStream, tmpTr);
  modelStream.clear();
  ReadParticles(partL, tmpTr);

  auto esti = std::make_shared<Estimator::SumMinLogLH>();
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
      partL, sqrtS4230._kin, 123);

  modelStream << modelSqrtS4230;
  xml_parser::read_xml(modelStream, tmpTr);
  modelStream.clear();

  //  sqrtS4230._data =
  //      std::make_shared<ComPWA::DataReader::RootReader>("data4230.root",
  //      "data");
  sqrtS4230._data = std::make_shared<Data>();
  sqrtS4230._mcSample = std::make_shared<Data>();
  ComPWA::Tools::GeneratePhsp(100000, sqrtS4230._gen, sqrtS4230._mcSample);

  // Construct intensity class from model string
  sqrtS4230._amp = IncoherentIntensity::Factory(
      partL, sqrtS4230._kin, tmpTr.get_child("IncoherentIntensity"));
  sqrtS4230._amp->GetParameters(fitPar);

  // We need to call this after the construction of the amplitude since
  // the variables are calculated that are needed by the amplitude
  sqrtS4230._mcPoints = std::make_shared<std::vector<dataPoint>>(
      sqrtS4230._mcSample->GetDataPoints(sqrtS4230._kin));
  sqrtS4230._amp->SetPhspSample(sqrtS4230._mcPoints, sqrtS4230._mcPoints);

  ComPWA::Tools::Generate(sqrtS4230._nEvents, sqrtS4230._kin, sqrtS4230._gen,
                          sqrtS4230._amp, sqrtS4230._data, sqrtS4230._mcSample);

  sqrtS4230._minLH = std::make_shared<Estimator::MinLogLH>(
      sqrtS4230._kin, sqrtS4230._amp, sqrtS4230._data, sqrtS4230._mcSample,
      sqrtS4230._mcSample, 0, 0);

  sqrtS4230._minLH->UseFunctionTree(true);
  sqrtS4230._minLH->GetTree()->Head()->SetName("logLH_sqrtS4230");
  esti->AddLogLh(sqrtS4230._minLH);

  //---------------------------------------------------
  // sqrtS = 4260
  //---------------------------------------------------
  energyPar sqrtS4260;
  cmsP4 = FourMomentum(0, 0, 0, 4.260);
  sqrtS4260._nEvents = 1000;

  sqrtS4260._kin = std::make_shared<HelicityKinematics>(partL, initialState,
                                                        finalState, cmsP4);
  sqrtS4260._gen = std::make_shared<ComPWA::Tools::RootGenerator>(
      partL, sqrtS4260._kin, 456);

  modelStream << modelSqrtS4260;
  xml_parser::read_xml(modelStream, tmpTr);

  //  sqrtS4260._data =
  //      std::make_shared<ComPWA::DataReader::RootReader>("data4230.root",
  //      "data");
  sqrtS4260._data = std::make_shared<Data>();
  sqrtS4260._mcSample = std::make_shared<Data>();
  ComPWA::Tools::GeneratePhsp(100000, sqrtS4260._gen, sqrtS4260._mcSample);

  // Construct intensity class from model string
  sqrtS4260._amp = IncoherentIntensity::Factory(
      partL, sqrtS4260._kin, tmpTr.get_child("IncoherentIntensity"));
  sqrtS4260._amp->GetParameters(fitPar);

  // We need to call this after the construction of the amplitude since
  // the variables are calculated that are needed by the amplitude
  sqrtS4260._mcPoints = std::make_shared<std::vector<dataPoint>>(
      sqrtS4260._mcSample->GetDataPoints(sqrtS4260._kin));
  sqrtS4260._amp->SetPhspSample(sqrtS4260._mcPoints, sqrtS4260._mcPoints);

  ComPWA::Tools::Generate(sqrtS4260._nEvents, sqrtS4260._kin, sqrtS4260._gen,
                          sqrtS4260._amp, sqrtS4260._data, sqrtS4260._mcSample);

  sqrtS4260._minLH = std::make_shared<Estimator::MinLogLH>(
      sqrtS4260._kin, sqrtS4260._amp, sqrtS4260._data, sqrtS4260._mcSample,
      sqrtS4260._mcSample, 0, 0);

  sqrtS4260._minLH->UseFunctionTree(true);
  sqrtS4260._minLH->GetTree()->Head()->SetName("logLH_sqrtS4260");
  esti->AddLogLh(sqrtS4260._minLH);

  //---------------------------------------------------
  // sqrtS = 4340
  //---------------------------------------------------

  //---------------------------------------------------
  // Run fit
  //---------------------------------------------------
  esti->UseFunctionTree(true);
  LOG(debug) << esti->GetTree()->Head()->Print(25);
  LOG(info) << "Fit parameter list: " << fitPar.to_str();
  auto minuitif = new Optimizer::Minuit2::MinuitIF(esti, fitPar);
  minuitif->SetHesse(true);
  auto result = minuitif->exec(fitPar);

  result->Print();

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
  sqrtS4230._pl = std::make_shared<RootPlot>(sqrtS4230._kin);
  sqrtS4230._pl->SetData(sqrtS4230._data);
  sqrtS4230._pl->SetPhspSample(sqrtS4230._mcSample);
  sqrtS4230._pl->SetFitAmp(sqrtS4230._amp);
  sqrtS4230._pl->AddComponent("f0(980)", "f0_980");
  sqrtS4230._pl->AddComponent("Zc(3900)_JpsiPiPlus + Zc(3900)_JpsiPiMinus",
                              "Zc3900");
  sqrtS4230._pl->Write("sqrtS4230", "plot.root", "RECREATE");

  sqrtS4260._pl = std::make_shared<RootPlot>(sqrtS4260._kin);
  sqrtS4260._pl->SetData(sqrtS4260._data);
  sqrtS4260._pl->SetPhspSample(sqrtS4260._mcSample);
  sqrtS4260._pl->SetFitAmp(sqrtS4260._amp);
  sqrtS4260._pl->AddComponent("f0(980)", "f0_980");
  sqrtS4260._pl->Write("sqrtS4260", "plot.root", "UPDATE");

  LOG(info) << "Done";

  return 0;
}
