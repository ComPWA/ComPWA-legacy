/*! test program : R -> XYZ
* @file testPWA.cpp
* Fit application for R -> \sum_i(R_i + X/Y/Z) -> X + Y + Z
*/

// Standard header files go here
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

// Root header files go here
#include "TFile.h"

// Boost header files
#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/xml_parser.hpp>

// Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/TableFormater.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"

#include "Physics/ParticleList.hpp"
#include "Physics/HelicityFormalism.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/CoherentIntensity.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"

#include "Tools/RootGenerator.hpp"
#include "Tools/Generate.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/DalitzPlot.hpp"
#include "Tools/ParameterTools.hpp"

using namespace std;
using namespace ComPWA;
using namespace ComPWA::DataReader;
using namespace ComPWA::Physics;
using namespace ComPWA::Physics::HelicityFormalism;

// Enable serialization of MinuitResult. For some reason has to be outside
// any namespaces.
BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)

/*****************************************************************************
 *
 * The main function.
 *
 *****************************************************************************/

int main(int argc, char **argv) {

  std::cout << " hello world" << std::endl;
  // read configuration for program setting/arguments
  namespace po = boost::program_options;

  std::string config_file;
  po::options_description generic("Generic options");
  generic.add_options()("help,h", "produce help message");
  generic.add_options()(
      "config,c", po::value<string>(&config_file)->default_value("config.cfg"),
      "name of a file of a configuration.");

  int seed;
  int numEvents;
  std::string logLevel;
  bool disableFileLog;
  std::string outputDir;
  std::string outputFileName;
  po::options_description config("General settings");
  config.add_options()(
      "logLevel", po::value<std::string>(&logLevel)->default_value("trace"),
      "set log level: error|warning|info|debug");
  config.add_options()("disableFileLog",
                       po::value<bool>(&disableFileLog)->default_value(0),
                       "write log file? 0/1");
  config.add_options()("nEvents",
                       po::value<int>(&numEvents)->default_value(100),
                       "set of events per fit");
  config.add_options()(
      "seed", po::value<int>(&seed)->default_value(110),
      "set random number seed; default is to used unique seed for every job");
  config.add_options()(
      "outputFile",
      po::value<std::string>(&outputFileName)->default_value("out"),
      "set output file name x. The files x.root and .log are created");
  config.add_options()("outputDir",
                       po::value<std::string>(&outputDir)->default_value("./"),
                       "set output directory");

  std::string dataFile, dataFileTreeName;
  std::string mcFile, mcFileTreeName;
  //  std::string bkgFile, bkgFileTreeName;
  std::string fitModelFile;
  std::string inputResult;
  po::options_description config_input("Input/Generate settings");
  config_input.add_options()(
      "dataFile", po::value<std::string>(&dataFile)->default_value(""),
      "set input data; leave blank for toy mc generation");
  config_input.add_options()(
      "dataFileTreeName",
      po::value<std::string>(&dataFileTreeName)->default_value("tree"),
      "set tree name in data file");
  config_input.add_options()("mcFile",
                             po::value<std::string>(&mcFile)->default_value(""),
                             "set input MC(phsp)");
  config_input.add_options()(
      "mcFileTreeName",
      po::value<std::string>(&mcFileTreeName)->default_value("tree"),
      "set tree name in mc file");
  //  config_input.add_options()(
  //      "bkgFile", po::value<std::string>(&bkgFile)->default_value("");
  //      "set bkg sample; leave blank if we fit data");
  //  config_input.add_options()(
  //      "bkgFileTreeName",
  //      po::value<std::string>(&bkgFileTreeName)->default_value("tree"),
  //      "set tree name in bkg file");
  config_input.add_options()(
      "fitModelFile",
      po::value<std::string>(&fitModelFile)->default_value("testModel.xml"),
      "set XML model file of data distribution");
  config_input.add_options()(
      "inputResult", po::value<std::string>(&inputResult)->default_value(""),
      "input file for MinuitResult");

  unsigned int mcPrecision;
  unsigned int ampMcPrecision;
  bool useMinos;
  bool useHesse;
  bool useRandomStartValues;
  //  bool calculateInterference;
  std::string fitMethod;
  std::string ampOption;
  po::options_description config_fit("Fit Settings");
  config_fit.add_options()(
      "mcPrecision",
      po::value<unsigned int>(&mcPrecision)->default_value(10000),
      "Precision form MC integration/cross section normalization. Number of "
      "events of phase space MC sample.");
  config_fit.add_options()(
      "ampMcPrecision",
      po::value<unsigned int>(&ampMcPrecision)->default_value(10000),
      "Precision for amplitude/breit-wigner normalization. Number of events of "
      "toy phase sample.");
  config_fit.add_options()(
      "fitMethod",
      po::value<std::string>(&fitMethod)->default_value("plotOnly"),
      "choose between 'tree', 'amplitude' and 'plotOnly'");
  config_fit.add_options()("useMinos",
                           po::value<bool>(&useMinos)->default_value(0),
                           "Run MINOS for each parameter");
  config_fit.add_options()("useHesse",
                           po::value<bool>(&useHesse)->default_value(1),
                           "Run HESSE after MIGRAD");
  config_fit.add_options()(
      "useRandomStartValues",
      po::value<bool>(&useRandomStartValues)->default_value(0),
      "Randomize start values");

  bool enablePlotting;
  bool DalitzPlotSample;
  int plotSize;
  int plotNBins;
  bool plotAmplitude;
  po::options_description config_plot("Plot Settings");
  config_plot.add_options()("plotting",
                            po::value<bool>(&enablePlotting)->default_value(1),
                            "Enable/Disable plotting");
  config_plot.add_options()(
      "DalitzPlot", po::value<bool>(&DalitzPlotSample)->default_value(1),
      "Eanble/Disable plotting of data");
  config_plot.add_options()("plotSize",
                            po::value<int>(&plotSize)->default_value(100000),
                            "Size of sample for amplitude plot");
  config_plot.add_options()("plotNBins",
                            po::value<int>(&plotNBins)->default_value(100),
                            "Number of bins");
  config_plot.add_options()("plotAmplitude",
                            po::value<bool>(&plotAmplitude)->default_value(1),
                            "Enable/Disable plotting of amplitude");

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config).add(config_input);
  cmdline_options.add(config_fit).add(config_plot);
  po::options_description config_file_options;
  config_file_options.add(config).add(config_input);
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

  Logging log(logFileName, boost::log::trivial::info); // initialize logging
  log.setLogLevel(logLevel);
  LOG(info) << " testPWA program for R-> XYZ";
  LOG(info) << " create/modify config.cfg to change configuration";
  LOG(info) << " create/modify testModel.xml to change the pwa amplitude model";
  LOG(info) << " create/modify minuitstragety.xml to change the Minuit2 "
               "stragety in minimization";

  // check configuration
  assert(!outputDir.empty());
  if (fitModelFile.empty() && fitMethod != "plotOnly") {
    LOG(error) << "Fit model file not set";
  }
  if (fitMethod != "tree" && fitMethod != "amplitude" &&
      fitMethod != "plotOnly") {
    LOG(error) << "Unknown fit method: " << fitMethod;
    return 1;
  }
  po::notify(vm);

  boost::property_tree::ptree fitModelTree;
  boost::property_tree::xml_parser::read_xml(fitModelFile, fitModelTree);
  // fill PartList
  auto partList = std::make_shared<ComPWA::PartList>();
  ReadParticles(partList, fitModelTree);

  //---------------------------------------------------
  // 1) Create Kinematics object
  //---------------------------------------------------
  auto fitModelKin = std::make_shared<HelicityKinematics>(
      partList, fitModelTree.get_child("HelicityKinematics"));
  std::cout << " step 1) ok" << std::endl;

  //---------------------------------------------------
  // 2) Create Intensity Class
  // this must before any generation/data,mc import
  //---------------------------------------------------
  auto intens = std::make_shared<IncoherentIntensity>(
      partList, fitModelKin, fitModelTree.get_child("Intensity"));
  boost::property_tree::ptree intens_tree = intens->save();
  std::cout << " write intensity tree " << std::endl;
  boost::property_tree::xml_parser::write_xml(std::cout, intens_tree);
  std::cout << " write intensity ok " << std::endl;
  // LOG(info) << "Subsystems used by fitModel:";
  // for (auto i : fitModelKin->subSystems()) {
  //  // Have to add " " here (bug in boost 1.59)
  //  LOG(info) << " " << i;
  //}

  std::cout << " step 2) ok" << std::endl;

  //---------------------------------------------------
  // 3) Generate a large phase space sample
  //---------------------------------------------------
  // initialize random generator using trueModel parameters
  std::shared_ptr<Generator> gen = std::shared_ptr<Generator>(
      new Tools::RootGenerator(partList, fitModelKin->initialState(),
                               fitModelKin->finalState(), seed));
  std::cout << " create generator ok " << std::endl;

  // WARNING: although it seems unlikely happen,
  // but one must create Intensity otherwise the HelicityKinematics object has
  // no subsystems
  // and then do the executions which needs the HelicityKinematics
  // otherwise, it will be wrong (althought it seems generator
  // will not use the intensity)

  // EFFICIENCY
  auto unitEff = std::shared_ptr<Efficiency>(new UnitEfficiency());
  std::cout << " unit efficiency ok " << std::endl;

  // Generation of toy phase space sample
  std::shared_ptr<Data> toySample(new RootReader()); // Toy sample
  ComPWA::Tools::generatePhsp(ampMcPrecision, gen, toySample);
  // set efficiency values for each event
  toySample->setEfficiency(fitModelKin, unitEff);
  auto toyPoints = std::make_shared<std::vector<DataPoint>>(
      toySample->dataPoints(fitModelKin));
  std::cout << " generate toySample ok " << std::endl;

  // Read Phase space sample which is used in PWA fit for xs section
  // normalization
  std::shared_ptr<Data> phspSample(new RootReader()); // phsp sample
  if (mcFile == "") {
    std::cout << " bbbb " << std::endl;
    // here using generation of Phase space instead of reading into
    ComPWA::Tools::generatePhsp(mcPrecision, gen, phspSample);
    std::cout << " generate phsp sample ok" << std::endl;
  } else {
    // read phase space sample
    std::shared_ptr<Data> phspSample(new RootReader(mcFile, mcFileTreeName));
    phspSample->reduceToPhsp(fitModelKin);
    std::cout << " import phspSample ok" << std::endl;
  }
  phspSample->reduceToPhsp(fitModelKin);
  phspSample->setEfficiency(fitModelKin, unitEff);
  auto phspPoints = std::make_shared<std::vector<DataPoint>>(
      phspSample->dataPoints(fitModelKin));
  std::cout << " phspsample ok" << std::endl;

  //
  intens->setPhspSample(phspPoints, toyPoints); // if only fit, this Set is not
                                                // necessary. If calculate
                                                // fitfraction, this set is
                                                // necessary
  std::cout << " set phspSample and toySample to intnesity ok " << std::endl;

  //---------------------------------------------------
  // 4) Read data which will be fit
  //---------------------------------------------------
  std::shared_ptr<Data> sample(new Data());
  if (dataFile != "") {
    std::shared_ptr<Data> inputData;
    //    int numSignalEvents =
    //        numEvents; // 445 data and 10000 mc for wrong matched signal
    LOG(info) << "Reading data file...";
    std::shared_ptr<Data> inD(new RootReader(dataFile, dataFileTreeName));
    inD->reduceToPhsp(fitModelKin);
    inD->setEfficiency(fitModelKin, unitEff);
    sample->append(*inD);
  } else {
    // generate a MC sample as data for test
    // int phspSize = numEvents/0.01;
    // std::shared_ptr<Data> phsp(new RooReader());
    // ComPWA::Tools::generatePhsp(phspSize, gen, phsp);

    LOG(info) << "Generating Sample";
    ComPWA::Tools::generate(numEvents, fitModelKin, gen, intens, sample,
                            std::shared_ptr<Data>(), std::shared_ptr<Data>());
    LOG(info) << "Sample size: " << sample->numEvents();
    LOG(info) << "Subsystems used by model: ";
    for (auto i : fitModelKin->subSystems()) {
      LOG(info) << " " << i;
    }
    std::stringstream s;
    s << "Printing the first 10 events of data sample:\n";
    for (int i = 0; (i < sample->numEvents() && i < 10); ++i) {
      DataPoint p;
      fitModelKin->convert(sample->event(i), p);
      s << p << "\n";
    }
    LOG(info) << s.str();

    std::cout << " bbb8 " << std::endl;
    sample->reduceToPhsp(fitModelKin);
    std::cout << " bbb9 " << std::endl;
    sample->setEfficiency(fitModelKin, unitEff);
  }

  //---------------------------------------------------
  // 5) Fit the model to the data and print the result
  //---------------------------------------------------

  LOG(info) << "==== FIT";
  LOG(info) << "MC precision: " << mcPrecision;
  LOG(info) << "AmpMc precision: " << ampMcPrecision;
  LOG(info) << "Fit method: " << fitMethod;
  LOG(info) << "Fit model file: " << fitModelFile;
  LOG(info) << "Use Hesse: " << useHesse;
  LOG(info) << "Use MINOS: " << useMinos;

  LOG(info) << "==== PLOTTING";
  LOG(info) << "eanble plotting: " << enablePlotting;
  if (enablePlotting) {
    LOG(info) << "plotting size: " << plotSize;
    LOG(info) << "number of bins: " << plotNBins;
    LOG(info) << "DalitzPlotSample: " << DalitzPlotSample;
    LOG(info) << "plotAmplitude: " << plotAmplitude;
  }
  LOG(info) << "=========================================";

  ParameterList fitPar;
  intens->parameters(fitPar);
  LOG(debug) << "Fit parameters: " << fitPar.to_str();
  // Set start error of 0.05 for parameters, run Minos?

  //=== Constructing likelihood
  auto esti = std::make_shared<Estimator::MinLogLH>(
      fitModelKin, intens, sample, toySample, phspSample, 0, 0);

  if (fitMethod == "tree") {
    esti->UseFunctionTree(true);
    esti->tree()->parameter();
    LOG(debug) << esti->tree()->head()->print(25);
  }

  if (useRandomStartValues)
    randomStartValues(fitPar);
  LOG(debug) << "Initial LH = " << esti->controlParameter(fitPar) << ".";

  // Set start error of 0.05 for parameters
  setErrorOnParameterList(fitPar, 0.05, useMinos);
  std::cout << " setErrorOnParameterList ok" << std::endl;

  auto minuitif = new Optimizer::Minuit2::MinuitIF(esti, fitPar);
  minuitif->setUseHesse(useHesse);
  std::cout << " setuseHesse ok " << std::endl;

  // STARTING MINIMIZATION
  std::shared_ptr<FitResult> result;
  std::cout << " begin exec " << std::endl;
  result = minuitif->exec(fitPar);
  std::cout << " minuitif->exec ok " << std::endl;

  // FIT RESULT
  ParameterList finalParList;
  //  std::cout << " begin get finalParameters " << std::endl;
  //  finalParList = result->finalParameters();
  //  std::cout << " fit ok" << std::endl;
  return 0;

  //  //calculation of fit fractions
  //  std::vector<std::pair<std::string, std::string>> fitComponents;
  //  fitComponents.push_back(
  //      std::pair<std::string, std::string>("phi(1020)", "D0toKSK+K-"));
  //  fitComponents.push_back(
  //      std::pair<std::string, std::string>("a0(980)0", "D0toKSK+K-"));
  //  fitComponents.push_back(
  //      std::pair<std::string, std::string>("a0(980)+", "D0toKSK+K-"));
  //  fitComponents.push_back(
  //      std::pair<std::string, std::string>("a2(1320)-", "D0toKSK+K-"));
  //  ParameterList ff = Tools::CalculateFitFractions(
  //      fitModelKin, intens->component("D0toKSK+K-"), toyPoints,
  //      fitComponents);
  //  result->setFitFractions(ff);
  //

  result->print();

  //---------------------------------------------------
  // 5.1) Save the fit result
  //---------------------------------------------------
  std::ofstream ofs(fileNamePrefix + std::string("-fitResult.xml"));
  boost::archive::xml_oarchive oa(ofs);
  oa << BOOST_SERIALIZATION_NVP(result);

  if (!disableFileLog) { // write fit result
    // Save final amplitude
    boost::property_tree::ptree ptout;
    ptout.add_child("IncoherentIntensity", intens->save());
    boost::property_tree::xml_parser::write_xml(
        fileNamePrefix + std::string("-Model.xml"), ptout, std::locale());
  } else if (!inputResult.empty()) {
    LOG(info) << "Reading MinuitResult from " << inputResult;
    std::shared_ptr<Optimizer::Minuit2::MinuitResult> inResult;
    std::ifstream ifs(inputResult);
    if (!ifs.good()) {
      throw std::runtime_error("input stream not good!");
    }
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(inResult);
    ifs.close();
    // if (calculateInterference) inResult->setCalcInterference(0);
    // inResult->SetFitFractions(ParameterList()); //reset fractions list
    result = inResult;
    result->print();
  }

  // Fill final parameters if minimization was not run
  if (!finalParList.numParameters()) {
    ParameterList tmpList;
    // Amplitude::FillAmpParameterToList(ampVec, tmpList);
    finalParList.DeepCopy(tmpList);
  }

  //---------------------------------------------------
  // 6) Plot data sample and intensity
  //---------------------------------------------------
  if (enablePlotting) {
    ComPWA::Tools::DalitzPlot pl(fitModelKin, fileNamePrefix, plotNBins);
    // set data sample
    pl.setData(sample);
    // set phsp sample
    pl.setPhspData(phspSample);
    // set amplitude
    pl.setFitAmp(intens, "", kBlue - 4);
    // select components to plot
    //    pl.drawComponent("D0toKsK+K-", "Signal", kGreen);

    // Fill histograms and create canvas
    pl.plot();
  }
  LOG(info) << "FINISHED";

  //// Exit code is exit code of fit routine. 0 is good/ 1 is bad
  if (result)
    return result->hasFailed();

  return 0;
}
