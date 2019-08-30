// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Executable for Dalitz plot analysis of the decay D0 -> K_S K+ K-
///

// Standard header files go here
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/export.hpp>

// Root header files go here
#include "TFile.h"

// Core header files go here
#include "Core/Event.hpp"
#include "Core/FitParameter.hpp"
#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/Logging.hpp"
#include "Core/Particle.hpp"
#include "Core/TableFormatter.hpp"

#include "Data/CorrectionTable.hpp"
#include "Data/DataCorrection.hpp"
#include "Data/Generate.hpp"
#include "Data/Root/RootDataIO.hpp"
#include "Data/Root/RootEfficiency.hpp"
#include "Data/Root/RootGenerator.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"

// HelicityFormlism
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Tools/FitFractions.hpp"

#include "Tools/Plotting/DalitzPlot.hpp"
#include "systematics.hpp"

using namespace std;
using namespace ComPWA;
using ComPWA::Data::DataSet;
using ComPWA::Data::MomentumCorrection;
using ComPWA::Data::Root::RootDataIO;
using ComPWA::Data::Root::RootEfficiency;
using namespace ComPWA::Physics;
using namespace ComPWA::Physics::HelicityFormalism;

// Enable serialization of MinuitResult. For some reason has to be outside
// any namespaces.
BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)

void applyCorrection(std::vector<Event> &sample,
                     MomentumCorrection &correction) {
  for (auto &x : sample) {
    x.Weight *= correction.correction(x);
  }
}

///
/// Dalitz plot analysis of the decay D0->K_S0 K_ K-.
///
int main(int argc, char **argv) {

  // Read command line input
  namespace po = boost::program_options;

  std::string config_file;

  po::options_description generic("Generic options");

  generic.add_options()("help,h", "produce help message");
  generic.add_options()(
      "config,c", po::value<string>(&config_file)->default_value("config.cfg"),
      "name of a file of a configuration.");

  int seed;      // initial seed
  int numEvents; // data size to be generated
  bool smearNEvents;
  std::string logLevel;
  bool disableFileLog;
  std::string outputDir;
  std::string outputFileName;

  po::options_description config("General settings");

  config.add_options()("nEvents",
                       po::value<int>(&numEvents)->default_value(1000),
                       "set of events per fit");
  config.add_options()(
      "seed", po::value<int>(&seed)->default_value(12345),
      "set random number seed; default is to used unique seed for every job");
  config.add_options()(
      "logLevel", po::value<std::string>(&logLevel)->default_value("trace"),
      "set log level: error|warning|info|debug");
  config.add_options()("disableFileLog",
                       po::value<bool>(&disableFileLog)->default_value(0),
                       "write log file? 0/1");
  config.add_options()(
      "outputFile",
      po::value<std::string>(&outputFileName)->default_value("out"),
      "set output file name x. The files x.root and .log are created");
  config.add_options()("outputDir",
                       po::value<std::string>(&outputDir)->default_value("./"),
                       "set output directory");
  config.add_options()("smearNEvents",
                       po::value<bool>(&smearNEvents)->default_value(0),
                       "smear NEvents for sqrt(NEvents)");

  std::string dataFile, dataFileTreeName;
  std::string bkgFile, bkgFileTreeName;
  std::string trueModelFile, trueBkgModelFile;
  double trueSignalFraction;
  std::string trueAmpOption;
  bool resetWeights;
  std::string inputResult;

  po::options_description config_inout("Input/Generate settings");

  config_inout.add_options()(
      "dataFile", po::value<std::string>(&dataFile)->default_value(""),
      "set input data; leave blank for toy mc generation");
  config_inout.add_options()(
      "dataFileTreeName",
      po::value<std::string>(&dataFileTreeName)->default_value("data"),
      "set tree name in data file");
  config_inout.add_options()(
      "bkgFile", po::value<std::string>(&bkgFile)->default_value(""),
      "set bkg sample; leave blank if we fit data");
  config_inout.add_options()(
      "bkgFileTreeName",
      po::value<std::string>(&bkgFileTreeName)->default_value("data"),
      "set tree name in bkg file");
  config_inout.add_options()(
      "trueModelFile",
      po::value<std::string>(&trueModelFile)->default_value("model.xml"),
      "set XML model file of true distribution");
  config_inout.add_options()(
      "trueBkgFile",
      po::value<std::string>(&trueBkgModelFile)->default_value(""),
      "set XML model file of true distribution");
  config_inout.add_options()(
      "trueSignalFraction",
      po::value<double>(&trueSignalFraction)->default_value(1.),
      "add this number of background events: Nbkg=(1-f)*NSignal; only makes "
      "sense together with bkgFile; ");
  config_inout.add_options()(
      "trueAmpOption",
      po::value<std::string>(&trueAmpOption)->default_value(""),
      "Use nocorrection=2/tagged=1/untagged=0");
  config_inout.add_options()("resetWeights",
                             po::value<bool>(&resetWeights)->default_value(0),
                             "should we reset all weights to one?");
  config_inout.add_options()(
      "inputResult", po::value<std::string>(&inputResult)->default_value(""),
      "input file for MinuitResult");

  std::string efficiencyFile;
  std::string efficiencyObject;
  std::string phspEfficiencyFile;
  std::string phspEfficiencyFileTreeName;
  std::string phspEfficiencyFileTrueTreeName;
  bool applySysCorrection;
  po::options_description config_eff("Efficiency settings");
  config_eff.add_options()(
      "efficiencyFile",
      po::value<std::string>(&efficiencyFile)->default_value(""),
      "root file with binned efficiency class. If none is specified uniform "
      "efficiency is assumed.")("efficiencyObject",
                                po::value<std::string>(&efficiencyObject)
                                    ->default_value("efficiency_m23cosTheta23"),
                                "name of TObject inside efficiencyFile")(
      "phspEfficiencyFile",
      po::value<std::string>(&phspEfficiencyFile)->default_value(""),
      "root file with reconstructed unbinned data from flat PHSP sample. We "
      "use this for "
      "an unbinned efficiency correction. Leave blank to use binned "
      "efficiency")(
      "phspEfficiencyFileTreeName",
      po::value<std::string>(&phspEfficiencyFileTreeName)->default_value("mc"),
      "set tree name for unbinned efficiency correction")(
      "phspEfficiencyFileTrueTreeName",
      po::value<std::string>(&phspEfficiencyFileTrueTreeName)
          ->default_value(""),
      "set tree name for tree with true values")(
      "applySysCorrection",
      po::value<bool>(&applySysCorrection)->default_value(0),
      "apply correction for data/MC difference");

  unsigned int mcPrecision;
  unsigned int ampMcPrecision;
  bool usePreFitter;
  bool useMinos;
  bool useHesse;
  bool useRandomStartValues;
  bool calculateInterference;
  int fitFractionError;
  double fitSignalFraction;
  std::string fitModelFile;
  std::string fitBkgFile;
  std::string fittingMethod;
  std::string ampOption;
  double penaltyScale;
  po::options_description config_fit("Fit settings");
  config_fit.add_options()(
      "mcPrecision",
      po::value<unsigned int>(&mcPrecision)->default_value(100000),
      "Precision for MC integration and normalization")(
      "ampMcPrecision",
      po::value<unsigned int>(&ampMcPrecision)->default_value(0),
      "Precision for MC integration and normalization")(
      "fittingMethod",
      po::value<std::string>(&fittingMethod)->default_value("tree"),
      "choose between 'tree' and 'plotOnly'")(
      "useMinos", po::value<bool>(&useMinos)->default_value(0),
      "Run MINOS for each parameter")(
      "useHesse", po::value<bool>(&useHesse)->default_value(1),
      "Run HESSE after MIGRAD")(
      "fitModelFile",
      po::value<std::string>(&fitModelFile)->default_value("model.xml"),
      "Set XML model file of fit distribution")(
      "fitBkgFile", po::value<std::string>(&fitBkgFile)->default_value(""),
      "Set XML background model file")(
      "ampOption", po::value<std::string>(&ampOption)->default_value("none"),
      "Use none/tagged/taggedSym/untagged")(
      "fitSignalFraction",
      po::value<double>(&fitSignalFraction)->default_value(1.),
      "Signal fraction")("usePreFitter",
                         po::value<bool>(&usePreFitter)->default_value(0),
                         "Run Geneva as prefitter")(
      "useRandomStartValues",
      po::value<bool>(&useRandomStartValues)->default_value(0),
      "Randomize start values")(
      "penaltyScale", po::value<double>(&penaltyScale)->default_value(0))(
      "calculateInterference",
      po::value<bool>(&calculateInterference)->default_value(0),
      "Calculate interference terms")(
      "fitFractionError", po::value<int>(&fitFractionError)->default_value(0),
      "Fit fraction errors are calculated by Monte-Carlo approach. 0 = "
      "disabled, >0 = MC precision");

  bool enablePlotting;
  bool DalitzPlotSample;
  int plotSize;
  bool plotCorrectEfficiency;
  int plotNBins;
  bool plotAmplitude;
  po::options_description config_plot("Plot settings");
  config_plot.add_options()("plotting",
                            po::value<bool>(&enablePlotting)->default_value(1),
                            "Enable/Disable plotting")(
      "DalitzPlot", po::value<bool>(&DalitzPlotSample)->default_value(1),
      "Enable/Disable plotting of data")(
      "plotSize", po::value<int>(&plotSize)->default_value(100000),
      "Size of sample for amplitude plot")(
      "plotNBins", po::value<int>(&plotNBins)->default_value(100),
      "Number of bins")(
      "plotCorrectEfficiency",
      po::value<bool>(&plotCorrectEfficiency)->default_value(0),
      "Number of bins")("plotAmplitude",
                        po::value<bool>(&plotAmplitude)->default_value(1),
                        "Enable/Disable plotting of amplitude");

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config).add(config_inout);
  cmdline_options.add(config_eff).add(config_fit).add(config_plot);
  po::options_description config_file_options;
  config_file_options.add(config).add(config_inout).add(config_eff);
  config_file_options.add(config_fit).add(config_plot);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);

  po::positional_options_description p;
  store(po::command_line_parser(argc, argv)
            .options(cmdline_options)
            .positional(p)
            .run(),
        vm);
  notify(vm);
  if (vm.count("help")) {
    cout << cmdline_options << "\n";
    return 1;
  }

  ifstream ifs(config_file.c_str());
  if (ifs) {
    store(parse_config_file(ifs, config_file_options), vm);
    notify(vm);
  }
  std::string fileNamePrefix = outputDir + std::string("/") + outputFileName;
  std::string logFileName = fileNamePrefix + std::string(".log");
  if (disableFileLog)
    logFileName = "";

  Logging log(logFileName, logLevel); // initialize logging

  // check configuration
  assert(!outputDir.empty());
  if (trueModelFile.empty()) {
    LOG(ERROR) << "True model file not set";
    trueModelFile = fitModelFile;
  }
  if (fitModelFile.empty() && fittingMethod != "plotOnly") {
    LOG(ERROR) << "Fit model file not set";
  }
  if (fittingMethod != "tree" && fittingMethod != "amplitude" &&
      fittingMethod != "plotOnly") {
    LOG(ERROR) << "Unknown fitting method: " << fittingMethod;
    return 1;
  }
  po::notify(vm);

  if (ampMcPrecision == 0)
    ampMcPrecision = mcPrecision;

  boost::property_tree::ptree fitModelTree;
  boost::property_tree::xml_parser::read_xml(fitModelFile, fitModelTree);
  boost::property_tree::ptree trueModelTree;
  boost::property_tree::xml_parser::read_xml(trueModelFile, trueModelTree);

  auto fitModelPartL = std::make_shared<ComPWA::PartList>();
  auto trueModelPartL = std::make_shared<ComPWA::PartList>();
  ReadParticles(fitModelPartL, fitModelTree);
  ReadParticles(trueModelPartL, trueModelTree);

  auto Builder = ComPWA::Physics::IntensityBuilderXML();
  // initialize kinematics of decay
  auto fitModelKin = Builder.createHelicityKinematics(
      fitModelPartL, fitModelTree.get_child("HelicityKinematics"));
  auto trueModelKin = Builder.createHelicityKinematics(
      trueModelPartL, trueModelTree.get_child("HelicityKinematics"));

  // Initialize random generator using trueModel parameters
  auto randGen = Data::Root::RootUniformRealGenerator(seed);

  auto gen = Data::Root::RootGenerator(
      trueModelKin.getParticleStateTransitionKinematicsInfo());

  // Binned efficiency
  if (!efficiencyFile.empty()) {
    TFile *tf = new TFile(TString(efficiencyFile));
    TH1 *h_passed = (TH1 *)tf->Get(TString(efficiencyObject) + "_passed");
    TH1 *h_total = (TH1 *)tf->Get(TString(efficiencyObject) + "_total");
    RootEfficiency rootEff(h_passed, h_total);
    tf->Close();
    auto eff = std::make_shared<RootEfficiency>(rootEff);
  }

  // Unbinned efficiency
  std::vector<Event> phspData, toyPhspData;

  if (!phspEfficiencyFile.empty()) {
    // sample with accepted phsp events
    RootDataIO RootReader(phspEfficiencyFileTreeName, mcPrecision);
    phspData = RootReader.readData(phspEfficiencyFile);
    phspData = Data::reduceToPhaseSpace(phspData, trueModelKin);
  }

  // Generation of toy phase space sample
  toyPhspData = ComPWA::Data::generatePhsp(mcPrecision, gen, randGen);
  if (phspData.size() == 0) {
    // in case no phase space sample is provided we use the toy sample
    phspData = toyPhspData;
  }

  Builder = Physics::IntensityBuilderXML(phspData);

  // Construct intensities
  auto intens = Builder.createIntensity(fitModelPartL, fitModelKin,
                                        fitModelTree.get_child("Intensity"));

  auto trueIntens = Builder.createIntensity(
      trueModelPartL, trueModelKin, trueModelTree.get_child("Intensity"));

  //======================= READING DATA =============================
  // Sample is used for minimization
  std::vector<Event> sample;

  // Temporary samples for signal and background
  std::vector<Event> inputData, inputBkg;

  // Smear total number of events
  // if (smearNEvents && numEvents > 0)  // I don't know what this is for...
  //  numEvents += (int)gen->gauss(0, sqrt(numEvents));

  if (!dataFile.empty()) {
    int numSignalEvents = numEvents;
    //    if (inputBkg)
    //      numSignalEvents -= inputBkg->getNEvents();
    LOG(INFO) << "Reading data file...";
    RootDataIO RR(dataFileTreeName);
    auto inD = RR.readData(dataFile);
    inD = Data::reduceToPhaseSpace(inD, trueModelKin);
    if (resetWeights) {
      for (auto &x : inD) {
        x.Weight = 1.0;
      }
    }
    // I dont understand the reason behind this...
    // This functionality was also removed. It can be put back in, if its really
    // needed
    // inputData = Data::randomSubSet(inputData, trueModelKin, numSignalEvents,
    // randGen);
    inputData = inD;

    sample = inputData;
  }

  //========== Generation of data sample ===========
  // generation with unbinned efficiency correction - use full sample
  if (!phspEfficiencyFile.empty() &&
      (inputData.size() == 0 ||
       inputBkg.size() == 0)) { // only if inputData is empty, right??
    // assume that efficiency of hit&miss is large 0.5%
    double phspSampleSize = numEvents / 0.005;

    RootDataIO RootReader(phspEfficiencyFileTreeName, phspSampleSize);
    auto fullPhsp = RootReader.readData(phspEfficiencyFile);
    if (applySysCorrection) {
      MomentumCorrection *trkSys =
          ComPWA::Data::getTrackingCorrection(fitModelPartL);
      MomentumCorrection *pidSys =
          ComPWA::Data::getPidCorrection(fitModelPartL);
      trkSys->print();
      pidSys->print();
      applyCorrection(fullPhsp, *trkSys);
      applyCorrection(fullPhsp, *pidSys);
    }
    std::vector<Event> fullTruePhsp; // Not used...
    if (!phspEfficiencyFileTrueTreeName.empty()) {
      RootDataIO RootReader2(phspEfficiencyFileTrueTreeName, phspSampleSize);
      fullTruePhsp = RootReader2.readData(phspEfficiencyFile);
      // Data* fullPhspRed = fullPhsp->EmptyClone();
      // Data* fullPhspTrueRed = fullTruePhsp->EmptyClone();
      // rndReduceSet(phspSampleSize,gen,fullPhsp.get(),fullPhspRed,
      // fullTruePhsp.get(),fullPhspTrueRed);
      // fullPhsp = std::shared_ptr<Data>(fullPhspRed);
      // fullTruePhsp = std::shared_ptr<Data>(fullPhspTrueRed);
    } else {
      // Data* fullPhspRed = fullPhsp->EmptyClone();
      // rndReduceSet(phspSampleSize,gen,fullPhsp.get(),fullPhspRed);
      // fullPhsp = std::shared_ptr<Data>(fullPhspRed);
    }
  }
  if (inputData.size() == 0) {
    LOG(INFO) << "Generating sample!";
    //    ComPWA::Tools::generate(numEvents, trueModelKin, gen, trueIntens,
    //    sample, phspData, toyPhspData);
    // Pass Null pointers for phsp sample which are than generated on the fly
    sample = ComPWA::Data::generate(numEvents, trueModelKin, gen, trueIntens,
                                    randGen);
    LOG(INFO) << "Sample size: " << sample.size();
  }
  // Reset phsp sample to save memory
  // run.SetPhspSample(std::shared_ptr<Data>());

  LOG(INFO) << "Subsystems used by true model:";
  for (auto i : trueModelKin.subSystems()) {
    LOG(INFO) << i;
  }
  std::stringstream s;
  s << "Printing the first 10 events of data sample:\n";
  for (unsigned int i = 0; (i < sample.size() && i < 10); ++i) {
    DataPoint p = trueModelKin.convert(sample[i]);
    s << p << "\n";
  }
  LOG(INFO) << s.str();

  // sample->writeData("test.root","tr");
  sample = Data::reduceToPhaseSpace(sample, trueModelKin);

  LOG(INFO) << "================== SETTINGS =================== ";

  LOG(INFO) << "==== GENERAL";
  LOG(INFO) << "Number of events: " << numEvents;
  LOG(INFO) << "Initial seed: " << seed;
  LOG(INFO) << "Log level: " << logLevel;
  LOG(INFO) << "Output file: " << fileNamePrefix << "-* [*.root,*.pdf,*.tex]";
  LOG(INFO) << "==== INPUT/GENERATION";
  if (!dataFile.empty()) { // read in data
    LOG(INFO) << "Data file: " << dataFile << " [tree=" << dataFileTreeName
              << " ]";
  } else { // generate signal
    LOG(INFO) << "true signal fraction: " << trueSignalFraction;
    LOG(INFO) << "True model file: " << trueModelFile;
    LOG(INFO) << "True amplitude option: " << trueAmpOption;
  }
  if (trueSignalFraction != 1.0) {
    if (!bkgFile.empty()) { // read in background
      LOG(INFO) << "Background file: " << bkgFile
                << " [tree=" << bkgFileTreeName << " ]";
    } else { // generate background
      LOG(INFO) << "True background model file: " << trueBkgModelFile;
    }
  }
  LOG(INFO) << "Total events in input sample: " << sample.size();

  LOG(INFO) << "==== EFFICIENCY";
  if (!efficiencyFile.empty())
    LOG(INFO) << "Binned Efficiency file: " << efficiencyFile
              << " [obj=" << efficiencyObject << "]";
  if (!phspEfficiencyFile.empty()) {
    LOG(INFO) << "PHSP data input file: " << phspEfficiencyFile
              << " [tree=" << phspEfficiencyFileTreeName << "]";
    if (!phspEfficiencyFileTrueTreeName.empty())
      LOG(INFO) << "True tree name: " << phspEfficiencyFileTrueTreeName;
    if (applySysCorrection)
      LOG(INFO) << "Unbinned efficiency is corrected "
                   "for tracking systematics!";
    LOG(INFO) << "Unbinned efficiency correction is used!";
  }
  if (phspEfficiencyFile.empty() && !efficiencyFile.empty())
    LOG(INFO) << "Binned efficiency correction is used!";
  if (!phspEfficiencyFile.empty() && !efficiencyFile.empty())
    LOG(INFO) << "No efficiency correction!";

  LOG(INFO) << "==== FIT";
  LOG(INFO) << "MC precision: " << mcPrecision;
  LOG(INFO) << "Amplitude MC precision: " << ampMcPrecision;
  LOG(INFO) << "Fitting method: " << fittingMethod;
  LOG(INFO) << "Fit signal fraction: " << fitSignalFraction;
  LOG(INFO) << "Amplitude option: " << ampOption;
  LOG(INFO) << "Fit model file: " << fitModelFile;
  if (!fitBkgFile.empty())
    LOG(INFO) << "Fit bkg model file: " << fitBkgFile;
  LOG(INFO) << "Using Geneva as pre fitter: " << usePreFitter;
  LOG(INFO) << "Use MINOS: " << useMinos;
  LOG(INFO) << "Accurate errors on fit fractions: " << fitFractionError;
  LOG(INFO) << "Penalty scale: " << penaltyScale;

  LOG(INFO) << "==== PLOTTING";
  LOG(INFO) << "enable plotting: " << enablePlotting;
  if (enablePlotting) {
    LOG(INFO) << "plotting size: " << plotSize;
    LOG(INFO) << "number of bins: " << plotNBins;
    LOG(INFO) << "Correct samples for efficiency: " << plotCorrectEfficiency;
    LOG(INFO) << "DalitzPlotSample: " << DalitzPlotSample;
    LOG(INFO) << "plotAmplitude: " << plotAmplitude;
  }
  LOG(INFO) << "===============================================";

  if (fittingMethod != "plotOnly") {
    //========================FITTING =====================
    // ComPWA::FunctionTree::createFitParameterList(truePar);

    //    std::cout << "phi = "
    //    <<Tools::Integral(intens->component("phi(1020)"), toyPoints,
    //                                 fitModelKin->phspVolume())
    //              << std::endl;
    //
    //    std::cout << "a0 = " <<Tools::Integral(intens->component("a0(980)0"),
    //    toyPoints,
    //                                 fitModelKin->phspVolume())
    //              << std::endl;
    //    std::cout << "total =
    //    "<<Tools::Integral(intens->component("D0toKSK+K-"), toyPoints,
    //                                 fitModelKin->phspVolume())
    //              << std::endl;
    //    std::cout << "ff= " << Tools::CalculateFitFraction(
    //                     fitModelKin, intens, toyPoints,
    //                     std::make_pair("phi(1020)", "D0toKSK+K-"))
    //              << std::endl;
    //=== Constructing likelihood
    auto esti = ComPWA::Estimator::createMinLogLHFunctionTreeEstimator(
        intens, Data::convertEventsToDataSet(sample, fitModelKin));

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

    Optimizer::Minuit2::MinuitIF minuitif(useHesse, useMinos);
    // Start minimization
    auto result = minuitif.optimize(esti.first, esti.second);

    // Calculation of fit fractions FEATURE CURRENTLY NOT AVAILABLE
    /*std::vector<std::pair<std::string, std::string>> fitComponents;
    fitComponents.push_back(
        std::pair<std::string, std::string>("phi(1020)", "D0toKSK+K-"));
    fitComponents.push_back(
        std::pair<std::string, std::string>("a0(980)0", "D0toKSK+K-"));
    fitComponents.push_back(
        std::pair<std::string, std::string>("a0(980)+", "D0toKSK+K-"));
    fitComponents.push_back(
        std::pair<std::string, std::string>("a2(1320)-", "D0toKSK+K-"));
    ParameterList ff = Tools::CalculateFitFractions(
        fitModelKin, intens->component("D0toKSK+K-"), toyPoints, fitComponents);
    Tools::CalcFractionError(
        fitPar,
        std::dynamic_pointer_cast<ComPWA::Optimizer::Minuit2::MinuitResult>(
            result)
            ->covarianceMatrix(),
        ff, fitModelKin, intens->component("D0toKSK+K-"), toyPoints, 20,
        fitComponents);
    result->setFitFractions(ff);*/

    // Print fit result
    LOG(INFO) << result;

    // Save fit result
    std::ofstream ofs(fileNamePrefix + std::string("-fitResult.xml"));
    boost::archive::xml_oarchive oa(ofs);
    oa << BOOST_SERIALIZATION_NVP(result);

    if (!inputResult.empty()) {
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
    }

    //======================= PLOTTING =============================
    if (enablePlotting) {
      std::vector<Event> pl_phspSample;
      LOG(INFO) << "Plotting results...";
      if (!phspEfficiencyFile.empty()) { // unbinned plotting
                                         // sample with accepted phsp events
        Data::Root::RootDataIO RR(phspEfficiencyFileTreeName, plotSize);
        pl_phspSample = RR.readData(phspEfficiencyFile);

        std::vector<Event> plotTruePhsp; // Unused again...
        if (!phspEfficiencyFileTrueTreeName.empty()) {
          Data::Root::RootDataIO RR(phspEfficiencyFileTrueTreeName, plotSize);
          plotTruePhsp = RR.readData(phspEfficiencyFile);
        }
      } else {
        pl_phspSample = ComPWA::Data::generatePhsp(plotSize, gen, randGen);
      }
      // reduce sample to phsp
      pl_phspSample = Data::reduceToPhaseSpace(pl_phspSample, trueModelKin);
      // pl_phspSample->setEfficiency(trueModelKin, eff); efficiency should be
      // included already

      //-------- Instance of DalitzPlot
      ComPWA::Tools::Plotting::DalitzPlot pl(fitModelKin, fileNamePrefix,
                                             plotNBins);
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
}
