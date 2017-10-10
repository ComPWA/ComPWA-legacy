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
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/ParticleList.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/RunManager.hpp"
#include "Tools/DalitzPlot.hpp"
#include "Tools/RootPlot.hpp"
#include "Tools/ParameterTools.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Estimator/MinLogLH/SumMinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

using namespace boost::property_tree;

using namespace ComPWA;
using namespace ComPWA::Tools;
using namespace ComPWA::Physics::HelicityFormalism;
using ComPWA::Physics::HelicityFormalism::IncoherentIntensity;
using ComPWA::Optimizer::Minuit2::MinuitResult;

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
  Logging log("DalitzFit-log.txt", boost::log::trivial::debug);

  ptree tmpTr;

  // List with all particle information needed
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, defaultParticleList);
  xml_parser::read_xml("particles.xml", tmpTr);
  ReadParticles(partL, tmpTr);

  auto esti = std::make_shared<Estimator::SumMinLogLH>();
  ParameterList fitPar;

  //---------------------------------------------------
  // sqrtS = 4230
  //---------------------------------------------------
  energyPar sqrtS4230;
  sqrtS4230._nEvents = 1000;
  xml_parser::read_xml("model4230.xml", tmpTr);
  sqrtS4230._kin = std::make_shared<HelicityKinematics>(
      partL, tmpTr.get_child("HelicityKinematics"));
    sqrtS4230._gen =
      std::make_shared<ComPWA::Tools::RootGenerator>(partL, sqrtS4230._kin);
  
  xml_parser::read_xml("model4230.xml", tmpTr);

  //  sqrtS4230._data =
  //      std::make_shared<ComPWA::DataReader::RootReader>("data4230.root",
  //      "data");
  sqrtS4230._data = std::make_shared<Data>();
  sqrtS4230._mcSample = std::make_shared<Data>();
  ComPWA::RunManager::genPhsp(100000, sqrtS4230._gen,
                              sqrtS4230._mcSample);
  
  // Construct intensity class from model string
  sqrtS4230._amp = IncoherentIntensity::Factory(
      partL, sqrtS4230._kin, tmpTr.get_child("IncoherentIntensity"));
  sqrtS4230._amp->GetParameters(fitPar);

  // We need to call this after the construction of the amplitude since
  // the variables are calculated that are needed by the amplitude
  sqrtS4230._mcPoints = std::make_shared<std::vector<dataPoint>>(
      sqrtS4230._mcSample->GetDataPoints(sqrtS4230._kin));
  sqrtS4230._amp->SetPhspSample(sqrtS4230._mcPoints, sqrtS4230._mcPoints);

  ComPWA::RunManager::gen(sqrtS4230._nEvents, sqrtS4230._kin, sqrtS4230._gen,
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
  sqrtS4260._nEvents = 1000;
  xml_parser::read_xml("model4230.xml", tmpTr);
  sqrtS4260._kin = std::make_shared<HelicityKinematics>(
      partL, tmpTr.get_child("HelicityKinematics"));
    sqrtS4260._gen =
      std::make_shared<ComPWA::Tools::RootGenerator>(partL, sqrtS4260._kin);
  
  xml_parser::read_xml("model4230.xml", tmpTr);

  //  sqrtS4260._data =
  //      std::make_shared<ComPWA::DataReader::RootReader>("data4230.root",
  //      "data");
  sqrtS4260._data = std::make_shared<Data>();
  sqrtS4260._mcSample = std::make_shared<Data>();
  ComPWA::RunManager::genPhsp(100000, sqrtS4260._gen,
                              sqrtS4260._mcSample);
  
  // Construct intensity class from model string
  sqrtS4260._amp = IncoherentIntensity::Factory(
      partL, sqrtS4260._kin, tmpTr.get_child("IncoherentIntensity"));
  sqrtS4260._amp->GetParameters(fitPar);

  // We need to call this after the construction of the amplitude since
  // the variables are calculated that are needed by the amplitude
  sqrtS4260._mcPoints = std::make_shared<std::vector<dataPoint>>(
      sqrtS4260._mcSample->GetDataPoints(sqrtS4260._kin));
  sqrtS4260._amp->SetPhspSample(sqrtS4260._mcPoints, sqrtS4260._mcPoints);

  ComPWA::RunManager::gen(sqrtS4260._nEvents, sqrtS4260._kin, sqrtS4260._gen,
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
  LOG(info) << "Fit parameter list: " << fitPar;
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
  sqrtS4230._pl =
      std::make_shared<RootPlot>(sqrtS4230._kin);
        sqrtS4230._pl->SetData(sqrtS4230._data);
  sqrtS4230._pl->SetPhspData(sqrtS4230._mcSample);
  sqrtS4230._pl->SetFitAmp(sqrtS4230._amp);
  sqrtS4230._pl->Write("4230", "plot.root", "RECREATE");

  sqrtS4260._pl =
      std::make_shared<RootPlot>(sqrtS4260._kin);
        sqrtS4260._pl->SetData(sqrtS4260._data);
  sqrtS4260._pl->SetPhspData(sqrtS4260._mcSample);
  sqrtS4260._pl->SetFitAmp(sqrtS4260._amp);
  sqrtS4260._pl->Write("4260", "plot.root", "UPDATE");

  LOG(info) << "Done";

  return 0;
}
