//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding functionality to generate set of
//events
//-------------------------------------------------------------------------------

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <memory>
#include <ctime>
#include <numeric>
#include <time.h>

// Root header files go here
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TMath.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

//Core header files go here
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Event.hpp"
#include "Core/Generator.hpp"
#include "Core/ProgressBar.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"

// ComPWA header files go here
#include "DataReader/Data.hpp"
#include "Core/Estimator.hpp"
#include "Optimizer/Optimizer.hpp"
#include "DataReader/RootReader/RootReader.hpp"
//#include "DataReader/JakeReader/JakeReader.hpp"
//#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
//#include "Optimizer/Geneva/GenevaIF.hpp"

#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "Physics/ParticleList.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/HelicityFormalism/IncoherentIntensity.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/progress.hpp>

#include "Examples/PythonDalitzFit/PythonFit.hpp"

#include "Tools/HistTools.hpp"
#include "Tools/ParameterTools.hpp"
#include "Tools/DalitzPlot.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/RunManager.hpp"
#include "Tools/DalitzPlot.hpp"

using namespace boost::log;

using namespace std;
using namespace ComPWA;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;
using ComPWA::Physics::HelicityFormalism::IncoherentIntensity;
using ComPWA::Optimizer::Minuit2::MinuitResult;
using ComPWA::DataReader::RootReader;



PythonFit::PythonFit() : argc(0), argv{}{
}


PythonFit::~PythonFit() {
    //LOG(debug) << "~RunManager: Last seed: " << gen_->getSeed();
}

int PythonFit::StartFit() {

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
	  std::shared_ptr<ComPWA::DataReader::Data> phspSample(new ComPWA::DataReader::Data());
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
	  std::shared_ptr<ComPWA::DataReader::Data> sample(new ComPWA::DataReader::Data());
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
