// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Executable for Dalitz plot analysis of the decay D0 -> K_S K+ K-
///
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/export.hpp>

// ComPWA header files
#include "Core/Event.hpp"
#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Data/CorrectionTable.hpp"
#include "Data/DataCorrection.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootDataIO.hpp"
#include "Data/Root/RootEfficiency.hpp"
#include "Data/Root/RootGenerator.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Plotting/DalitzPlot.hpp"
#include "systematics.hpp"

using namespace ComPWA;
using ComPWA::Data::DataSet;
using ComPWA::Data::MomentumCorrection;
using ComPWA::Data::Root::RootDataIO;
using namespace ComPWA::Physics;
using namespace ComPWA::Physics::HelicityFormalism;

namespace po = boost::program_options;

// Enable serialization of MinuitResult. For some reason has to be outside
// any namespaces.
// BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)

///
/// Dalitz plot analysis of the decay D0->K_S0 K_ K-.
///
int main(int argc, char **argv) {

  // --------- Configuration options using boost::program_options --------
  std::string config_file;
  int seed; // initial seed
  std::string logLevel;
  std::string pathPrefix;
  std::string logFileName;

  po::options_description config("Settings");
  config.add_options()("help,h", "produce help message");
  config.add_options()(
      "config,c",
      po::value<std::string>(&config_file)->default_value("config.cfg"),
      "Configuration file name");
  config.add_options()(
      "seed", po::value<int>(&seed)->default_value(12345),
      "Set random number seed (0 to use cput for initial seed)");
  config.add_options()(
      "logLevel", po::value<std::string>(&logLevel)->default_value("debug"),
      "Set log level: error|warning|info|debug(default)|trace");
  config.add_options()("logFileName",
                       po::value<std::string>(&logFileName)->default_value(""),
                       "Log file name (empty to disable logging to file)");
  config.add_options()("pathPrefix",
                       po::value<std::string>(&pathPrefix)->default_value("./"),
                       "Set prefix of output files");

  std::string dataFile, dataFileTreeName;
  int numEvents; // data size to be generated
  bool smearNEvents;
  bool resetWeights;
  std::string trueModelFile;
  std::string inputResult;
  std::string phspSampleFile;
  std::string phspSampleFileTreeName;
  std::string phspSampleFileTrueTreeName;
  bool phspSample_applySysCorrection;

  config.add_options()(
      "dataFile", po::value<std::string>(&dataFile)->default_value(""),
      "Data sample (leave empty to generate MC according to true model)");
  config.add_options()(
      "dataFileTreeName",
      po::value<std::string>(&dataFileTreeName)->default_value("data"),
      "Name of TTree of data sample");
  config.add_options()("nEvents",
                       po::value<int>(&numEvents)->default_value(1000),
                       "set of events per fit");
  config.add_options()("smearNEvents",
                       po::value<bool>(&smearNEvents)->default_value(0),
                       "Vary number of events by poission statistics");
  config.add_options()("resetWeights",
                       po::value<bool>(&resetWeights)->default_value(0),
                       "Reset existing weights in data sample");
  config.add_options()(
      "trueModelFile",
      po::value<std::string>(&trueModelFile)->default_value("model.xml"),
      "File for true model");
  config.add_options()("inputResult",
                       po::value<std::string>(&inputResult)->default_value(""),
                       "Read fit result from file");
  config.add_options()(
      "phspSampleFile",
      po::value<std::string>(&phspSampleFile)->default_value(""), "");
  config.add_options()(
      "phspSampleFileTreeName",
      po::value<std::string>(&phspSampleFileTreeName)->default_value("mc"),
      "Name of TTree of phsp sample (reconstructed position)");
  config.add_options()(
      "phspSampleFileTrueTreeName",
      po::value<std::string>(&phspSampleFileTrueTreeName)->default_value(""),
      "Name of TTree of phsp sample (generated position)");
  config.add_options()(
      "applySysCorrection",
      po::value<bool>(&phspSample_applySysCorrection)->default_value(0),
      "Apply systematic tracking and PID corrections");

  size_t mcPrecision;
  bool useMinos;
  bool useRandomStartValues;
  bool calculateInterference;
  int fitFractionError;
  std::string fitModelFile;
  bool enableFit;
  std::string fitResultFile;
  bool enablePlot;
  int plotSize;
  std::string plotFileName;

  po::options_description config_fit("Fit/Plot settings");
  config_fit.add_options()(
      "mcPrecision", po::value<size_t>(&mcPrecision)->default_value(100000),
      "Precision for MC integration and normalization");
  config_fit.add_options()("enableFit",
                           po::value<bool>(&enableFit)->default_value(1),
                           "Enable minimization of likelihood");
  config_fit.add_options()("useMinos",
                           po::value<bool>(&useMinos)->default_value(0),
                           "Run MINOS for each parameter");
  config_fit.add_options()(
      "fitModelFile",
      po::value<std::string>(&fitModelFile)->default_value("model.xml"),
      "Set XML model file of fit distribution");
  config_fit.add_options()(
      "useRandomStartValues",
      po::value<bool>(&useRandomStartValues)->default_value(0),
      "Randomize start values");
  config_fit.add_options()(
      "calculateInterference",
      po::value<bool>(&calculateInterference)->default_value(0),
      "Calculate interference terms");
  config_fit.add_options()(
      "fitFractionError", po::value<int>(&fitFractionError)->default_value(0),
      "Fit fraction errors are calculated by Monte-Carlo approach. 0 = "
      "disabled, >0 = MC precision");
  config_fit.add_options()(
      "fitResultFile", po::value<std::string>(&fitResultFile)->default_value("fitResult.xml"),
      "File to save fit result");
  config_fit.add_options()("enablePlot",
                           po::value<bool>(&enablePlot)->default_value(1),
                           "Enable/Disable plotting");
  config_fit.add_options()("plotSize",
                           po::value<int>(&plotSize)->default_value(100000),
                           "Number of events used to evaluate model");
  config_fit.add_options()(
      "plotFileName",
      po::value<std::string>(&plotFileName)->default_value("plot"),
      "File name for plots (obmit the suffix *.root)");

  po::options_description programOpts;
  programOpts.add(config).add(config_fit);

  // Read command line arguments
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, programOpts), vm);

  po::positional_options_description p;
  store(po::command_line_parser(argc, argv)
            .options(programOpts)
            .positional(p)
            .run(),
        vm);
  notify(vm);
  if (vm.count("help")) {
    std::cout << programOpts << "\n";
    return 1;
  }
  // Read config file
  std::ifstream ifs(config_file.c_str());
  if (ifs) {
    store(parse_config_file(ifs, programOpts), vm);
    notify(vm);
  }
  // --------- Configuration options END --------

  pathPrefix += "/"; // make sure '/' is the last char

  if (!logFileName.empty())
    logFileName = pathPrefix + logFileName;

  Logging log(logLevel, logFileName); // initialize logging

  // check configuration
  if (trueModelFile.empty())
    trueModelFile = fitModelFile;

  auto Builder = ComPWA::Physics::IntensityBuilderXML();

  // Build the true model
  boost::property_tree::ptree trueModelTree;
  boost::property_tree::xml_parser::read_xml(trueModelFile, trueModelTree);
  auto trueParticleList = std::make_shared<ComPWA::PartList>();
  ReadParticles(trueParticleList, trueModelTree);

  auto trueKinematics = Builder.createHelicityKinematics(
      trueParticleList, trueModelTree.get_child("HelicityKinematics"));

  // Build the fit model
  boost::property_tree::ptree fitModelTree;
  boost::property_tree::xml_parser::read_xml(fitModelFile, fitModelTree);
  auto fitParticleList = std::make_shared<ComPWA::PartList>();
  ReadParticles(fitParticleList, fitModelTree);

  auto fitKinematics = Builder.createHelicityKinematics(
      fitParticleList, fitModelTree.get_child("HelicityKinematics"));


  auto randGen = StdUniformRealGenerator(seed);
  std::mt19937 randGen2(seed);

  // Initialize random generator using trueModel parameters
  auto gen = Data::Root::RootGenerator(
      trueKinematics.getParticleStateTransitionKinematicsInfo());

  // Two samples of phsp events can be specified:
  // 1) A sample which included the reconstruction efficiency (e.g. a sample
  //    of events which have passed detector simulation, reconstruction and
  //    selection. From this sample we need:
  //    - the reconstructed phsp position -> phspSample
  //    - the generated phsp position -> phspSampleTrue
  //    This sample is used to correct the likelihood for efficiency and to
  //    properly generate events including detector effects.
  // 2) A toy phsp sample which is generated on the fly -> phspSampleToy
  // Samples 1) and 2) can be independent of each other. In the most simple case
  // the toy sample is used for all of them.
  std::vector<Event> phspSample, phspSampleTrue, phspSampleToy;

  if (!phspSampleFile.empty()) {
    // sample with accepted phsp events
    RootDataIO RootReader(phspSampleFileTreeName, mcPrecision);
    phspSample = RootReader.readData(phspSampleFile);
    phspSample = Data::reduceToPhaseSpace(phspSample, trueKinematics);
    if (!phspSampleFileTrueTreeName.empty()) {
      RootDataIO RootReaderTrue(phspSampleFileTrueTreeName, mcPrecision);
      phspSampleTrue = RootReaderTrue.readData(phspSampleFile);
    } else {
      phspSampleTrue = phspSample;
    }
    if (phspSample_applySysCorrection) {
      MomentumCorrection trkSys =
          ComPWA::Data::getTrackingCorrection(fitParticleList);
      LOG(INFO) << trkSys;
      std::for_each(phspSample.begin(), phspSample.end(),
                    [&trkSys](Event &ev) { ev.Weight *= trkSys(ev); });

      MomentumCorrection pidSys =
          ComPWA::Data::getPidCorrection(fitParticleList);
      LOG(INFO) << pidSys;
      std::for_each(phspSample.begin(), phspSample.end(),
                    [&pidSys](Event &ev) { ev.Weight *= pidSys(ev); });
    }
  }

  // Generation of toy phase space sample
  phspSampleToy = ComPWA::Data::generatePhsp(mcPrecision, gen, randGen);
  if (phspSample.size() == 0) {
    // in case no phase space sample is provided we use the toy sample
    phspSample = phspSampleTrue = phspSampleToy;
  }

  // TODO: Builder needs two samples here (see Issue #213)
  //  Builder = Physics::IntensityBuilderXML(phspSample, phspSampleToy);
  Builder = Physics::IntensityBuilderXML(phspSample);

  auto trueIntens = Builder.createIntensity(
      trueParticleList, trueKinematics, trueModelTree.get_child("Intensity"));
  LOG(INFO) << "Subsystems used by true model:";
  for (auto i : trueKinematics.subSystems()) {
    LOG(INFO) << i;
  }
  //======================= Reading/Generating data sample =====================
  // Smear total number of events
  if (smearNEvents && numEvents > 0) {
    numEvents = std::poisson_distribution<int>(numEvents)(randGen2);
  }

  std::vector<Event> sample;
    if (!dataFile.empty()) {
    LOG(INFO) << "Reading data sample from file...";
    RootDataIO RR(dataFileTreeName,
                  numEvents); // numEvents = -1 -> read all data
    sample = RR.readData(dataFile);
    sample = Data::reduceToPhaseSpace(sample, trueKinematics);
    if (resetWeights) {
      for (auto &x : sample) {
        x.Weight = 1.0;
      }
    }

    // Sample size is larger than requested -> select a random subset
    if(sample.size() > numEvents){
      sample.resize(numEvents);
      std::shuffle(sample.begin(), sample.end(), randGen2);
    }
  } else {
    LOG(INFO) << "Generating data sample...";
    sample = ComPWA::Data::generate(numEvents, trueKinematics, randGen,
                                    trueIntens, phspSample, phspSampleTrue);
  }
  sample = Data::reduceToPhaseSpace(sample, trueKinematics);

  std::stringstream s;
  s << "Printing the first 10 events of data sample:\n";
  for (unsigned int i = 0; (i < sample.size() && i < 10); ++i)
    s << trueKinematics.convert(sample.at(i)) << "\n";
  LOG(INFO) << s.str();

  LOG(INFO) << "================== SETTINGS =================== ";
  LOG(INFO) << "==== General";
  LOG(INFO) << "Initial seed: " << seed;
  LOG(INFO) << "Path prefix: " << pathPrefix;
  LOG(INFO) << "Log file (level): " << logFileName << " (" << logLevel << ")";
  LOG(INFO) << "==== Input samples and true model";
  LOG(INFO) << "Number of events: " << numEvents;
  if (!dataFile.empty()) // read in data
    LOG(INFO) << "Data file: " << dataFile << " [tree=" << dataFileTreeName
              << " ]";
  else // generate signal
    LOG(INFO) << "True model file: " << trueModelFile;
  LOG(INFO) << "Total events in input sample: " << sample.size();
  if (!phspSampleFile.empty()) {
    LOG(INFO) << "PHSP data input file: " << phspSampleFile
              << " [tree=" << phspSampleFileTreeName << "]";
    if (!phspSampleFileTrueTreeName.empty())
      LOG(INFO) << "TTree name of generated phsp position: "
                << phspSampleFileTrueTreeName;
    LOG(INFO) << "Correct tracking/pid systematics: "
              << phspSample_applySysCorrection;
  }
  LOG(INFO) << "==== FIT";
  LOG(INFO) << "Enable fit: " << enableFit;
  LOG(INFO) << "MC precision: " << mcPrecision;
  LOG(INFO) << "Fit model file: " << fitModelFile;
  LOG(INFO) << "Use MINOS: " << useMinos;
  LOG(INFO) << "Accurate errors on fit fractions: " << fitFractionError;
  LOG(INFO) << "==== PLOTTING";
  LOG(INFO) << "Enable plotting: " << enablePlot;
  if (enablePlot)
    LOG(INFO) << "Number of events to plot: " << plotSize;
  LOG(INFO) << "===============================================";

  // Fit model
  auto fitIntens = Builder.createIntensity(fitParticleList, fitKinematics,
                                           fitModelTree.get_child("Intensity"));
  ComPWA::Optimizer::Minuit2::MinuitResult result;

  if (enableFit) {
    //========================FITTING =====================

    // Construct likelihood
    auto esti = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
        fitIntens, Data::convertEventsToDataSet(sample, fitKinematics));

    for (auto x : esti.second)
      LOG(DEBUG) << x;

    std::cout.setf(std::ios::unitbuf);
    LOG(DEBUG) << esti.first.print(25);

    if (useRandomStartValues) {
      for (auto &x : esti.second) {
        if (!x.IsFixed) {
          double RandomValue(randGen());
          std::pair<double, double> bounds(-999, -999);
          if (x.HasBounds) {
            bounds = x.Bounds;
          }
          x.Value = bounds.first + RandomValue * (bounds.second - bounds.first);
          LOG(DEBUG) << x;
        }
      }
    }

    Optimizer::Minuit2::MinuitIF minuitif(true, useMinos);

    // Start minimization
    result = minuitif.optimize(esti.first, esti.second);

    // Calculation of fit fractions FEATURE CURRENTLY NOT AVAILABLE
    // std::vector<std::pair<std::string, std::string>> fitComponents;
    //    fitComponents.push_back(
    //        std::pair<std::string, std::string>("phi(1020)", "D0toKSK+K-"));
    //    fitComponents.push_back(
    //        std::pair<std::string, std::string>("a0(980)0", "D0toKSK+K-"));
    //    fitComponents.push_back(
    //        std::pair<std::string, std::string>("a0(980)+", "D0toKSK+K-"));
    //    fitComponents.push_back(
    //        std::pair<std::string, std::string>("a2(1320)-", "D0toKSK+K-"));
    //    ParameterList ff = Tools::CalculateFitFractions(
    //        fitModelKin, intens->component("D0toKSK+K-"), toyPoints,
    //        fitComponents);
    //    Tools::CalcFractionError(
    //        fitPar,
    //        std::dynamic_pointer_cast<ComPWA::Optimizer::Minuit2::MinuitResult>(
    //            result)
    //            ->covarianceMatrix(),
    //        ff, fitModelKin, intens->component("D0toKSK+K-"), toyPoints, 20,
    //        fitComponents);
    //    result->setFitFractions(ff);

    // Print fit result
    LOG(INFO) << result;

    // Save fit result
    std::ofstream ofs(pathPrefix + fitResultFile);
    boost::archive::xml_oarchive oa(ofs);
    oa << BOOST_SERIALIZATION_NVP(result);

  } else if (!inputResult.empty()) {
    // If no fit is performed we can read in an existing fit result, update the
    // intensity and plot it in the next step.
    LOG(INFO) << "Reading MinuitResult from " << inputResult;
    Optimizer::Minuit2::MinuitResult inResult;
    std::ifstream ifs(inputResult);
    if (!ifs.good())
      throw std::runtime_error("input stream not good!");
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(inResult);
    ifs.close();

    result = inResult;
    LOG(INFO) << result;
    // TODO: Update intensity
  }

  //======================= PLOTTING =============================
  if (enablePlot) {
    std::vector<Event> pl_phspSample;
    LOG(INFO) << "Plotting results...";
    if (!phspSampleFile.empty()) {
      Data::Root::RootDataIO RR(phspSampleFileTreeName, plotSize);
      pl_phspSample = RR.readData(phspSampleFile);

      std::vector<Event> plotTruePhsp; // Unused again...
      if (!phspSampleFileTrueTreeName.empty()) {
        Data::Root::RootDataIO RR(phspSampleFileTrueTreeName, plotSize);
        plotTruePhsp = RR.readData(phspSampleFile);
      }
    } else {
      pl_phspSample = ComPWA::Data::generatePhsp(plotSize, gen, randGen);
    }
    // reduce sample to phsp
    pl_phspSample = Data::reduceToPhaseSpace(pl_phspSample, trueKinematics);
    // pl_phspSample->setEfficiency(trueModelKin, eff); efficiency should be
    // included already

    //-------- Instance of DalitzPlot
    ComPWA::Tools::Plotting::DalitzPlot pl(fitKinematics,
                                           pathPrefix + plotFileName, 100);
    // set data sample
    pl.fillData(sample);
    // set phsp sample
    // pl.fillPhaseSpaceData(pl_phspSample, intens);
    // set amplitude
    /*pl.setFitAmp(intens, "", kBlue - 4);
    // select components to plot
    pl.drawComponent("D0toKSK+K-", "D0toKSK+K-", "Signal", kGreen);
    pl.drawComponent("a0(980)0", "D0toKSK+K-", "a_{0}(980)^{0}", kMagenta);
    pl.drawComponent("a0(980)+", "D0toKSK+K-", "a_{0}(980)^{+}",
                     kMagenta + 2);
    pl.drawComponent("phi(1020)", "D0toKSK+K-", "#phi(1020)", kMagenta + 4);
    pl.drawComponent("BkgD0toKSK+K-", "BkgD0toKSK+K-", "Background", kRed);*/

    // Fill histograms and create canvases
    pl.plot();
  }
  LOG(INFO) << "FINISHED!";

  // Exit code is exit code of fit routine. 0 is good/ 1 is bad
  return result.IsValid;
}
