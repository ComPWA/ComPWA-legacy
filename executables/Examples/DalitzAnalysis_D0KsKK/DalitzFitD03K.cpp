/*! Dalitz plot analysis of the decay D0->K_S0 K_ K-
 * @file DalitzFitApp.cpp
 * Fit application for D0->K_S0 K_ K- analysis. A model with BW and flatte formalism is used to
 * describe the data.
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

// Root header files go here
#include "TFile.h"

//Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/RunManager.hpp"
#include "Core/Efficiency.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Amplitude.hpp"
#include "Core/TableFormater.hpp"
#include "Core/AbsParameter.hpp"
#include "Core/Logging.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "DataReader/CorrectionTable.hpp"
#include "DataReader/DataCorrection.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Physics/DPKinematics/RootEfficiency.cpp"
#include "Physics/DPKinematics/SimpleEfficiency.hpp"
#include "Physics/DPKinematics/SimpleAngleEfficiency.hpp"
#include "Physics/DPKinematics/RootGenerator.hpp"

#include "PlotData.hpp"
#include "Tools.hpp"

using namespace std;
using namespace ComPWA;

BOOST_CLASS_EXPORT(Optimizer::Minuit2::MinuitResult)

/*****************************************************************************
 *
 * The main function.
 *
 *****************************************************************************/

int main(int argc, char **argv){
    
    // Read command line input
    namespace po = boost::program_options;
    
    std::string config_file;
    po::options_description generic("Generic options");
    generic.add_options()
    ("help,h", "produce help message")
    ("config,c", po::value<string>(&config_file)->default_value(""),
     "name of a file of a configuration.");
    
    int seed; //initial seed
    int numEvents;//data size to be generated
    bool smearNEvents;
    std::string logLevel;
    bool disableFileLog;
    std::string outputDir;
    std::string outputFileName;
    po::options_description config("General settings");
    config.add_options()
    ("nEvents", po::value<int>(&numEvents)->default_value(1000), "set of events per fit")
    ("seed", po::value<int>(&seed)->default_value(-1), "set random number seed; default is to used unique seed for every job")
    ("logLevel", po::value<std::string>(&logLevel)->default_value("info"), "set log level: error|warning|info|debug")
    ("disableFileLog", po::value<bool>(&disableFileLog)->default_value(0), "write log file? 0/1")
    ("outputFile", po::value<std::string>(&outputFileName)->default_value("out"),
     "set output file name x. The files x.root and .log are created")
    ("outputDir", po::value<std::string>(&outputDir)->default_value("./"),
     "set output directory")
    ("smearNEvents", po::value<bool>(&smearNEvents)->default_value(0), "smear NEvents for sqrt(NEvents)");
    
    std::string dataFile, dataFileTreeName;
    std::string bkgFile, bkgFileTreeName;
    std::string trueModelFile, trueBkgModelFile;
    double trueSignalFraction;
    std::string trueAmpOption;
    bool resetWeights;
    std::string inputResult;
    po::options_description config_inout("Input/Generate settings");
    config_inout.add_options()
    ("dataFile", po::value<std::string>(&dataFile)->default_value(""),
     "set input data; leave blank for toy mc generation")
    ("dataFileTreeName", po::value<std::string>(&dataFileTreeName)->default_value("data"),
     "set tree name in data file")
    ("bkgFile", po::value<std::string>(&bkgFile)->default_value(""),
     "set bkg sample; leave blank if we fit data")
    ("bkgFileTreeName", po::value<std::string>(&bkgFileTreeName)->default_value("data"),
     "set tree name in bkg file")
    ("trueModelFile", po::value<std::string>(&trueModelFile)->default_value("model.xml"), "set XML model file of true distribution")
    ("trueBkgFile", po::value<std::string>(&trueBkgModelFile)->default_value(""), "set XML model file of true distribution")
    ("trueSignalFraction", po::value<double>(&trueSignalFraction)->default_value(1.),
     "add this number of background events: Nbkg=(1-f)*NSignal; only makes sense together with bkgFile; ")
    ("trueAmpOption", po::value<std::string>(&trueAmpOption)->default_value(""), "Use nocorrection=2/tagged=1/untagged=0")
    ("resetWeights", po::value<bool>(&resetWeights)->default_value(0), "should we reset all weights to one?")
    ("inputResult", po::value<std::string>(&inputResult)->default_value(""),
     "input file for MinuitResult");
    
    std::string efficiencyFile;
    std::string efficiencyObject;
    std::string phspEfficiencyFile;
    std::string phspEfficiencyFileTreeName;
    std::string phspEfficiencyFileTrueTreeName;
    bool applySysCorrection;
    po::options_description config_eff("Efficiency settings");
    config_eff.add_options()
    ("efficiencyFile", po::value<std::string>(&efficiencyFile)->default_value(""),
     "root file with binned efficiency class. If none is specified uniform efficiency is assumed.")
    ("efficiencyObject", po::value<std::string>(&efficiencyObject)->default_value("efficiency_m23cosTheta23"),
     "name of TObject inside efficiencyFile")
    ("phspEfficiencyFile", po::value<std::string>(&phspEfficiencyFile)->default_value(""),
     "root file with reconstructed unbinned data from flat PHSP sample. We use this for "
     "an unbinned efficiency correction. Leave blank to use binned efficiency")
    ("phspEfficiencyFileTreeName", po::value<std::string>(&phspEfficiencyFileTreeName)->default_value("mc"),
     "set tree name for unbinned efficiency correction")
    ("phspEfficiencyFileTrueTreeName", po::value<std::string>(&phspEfficiencyFileTrueTreeName)->default_value(""),
     "set tree name for tree with true values")
    ("applySysCorrection", po::value<bool>(&applySysCorrection)->default_value(0),
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
    config_fit.add_options()
    ("mcPrecision", po::value<unsigned int>(&mcPrecision)->default_value(100000), "Precision for MC integration and normalization")
    ("ampMcPrecision", po::value<unsigned int>(&ampMcPrecision)->default_value(0), "Precision for MC integration and normalization")
    ("fittingMethod", po::value<std::string>(&fittingMethod)->default_value("plotOnly"), "choose between 'tree', 'amplitude' and 'plotOnly'")
    ("useMinos", po::value<bool>(&useMinos)->default_value(0), "Run MINOS for each parameter")
    ("useHesse", po::value<bool>(&useHesse)->default_value(1), "Run HESSE after MIGRAD")
    ("fitModelFile", po::value<std::string>(&fitModelFile)->default_value("model.xml"), "Set XML model file of fit distribution")
    ("fitBkgFile", po::value<std::string>(&fitBkgFile)->default_value(""), "Set XML background model file")
    ("ampOption", po::value<std::string>(&ampOption)->default_value("none"), "Use none/tagged/taggedSym/untagged")
    ("fitSignalFraction", po::value<double>(&fitSignalFraction)->default_value(1.),	"Signal fraction")
    ("usePreFitter", po::value<bool>(&usePreFitter)->default_value(0), "Run Geneva as prefitter")
    ("useRandomStartValues", po::value<bool>(&useRandomStartValues)->default_value(0), "Randomize start values")
    ("penaltyScale", po::value<double>(&penaltyScale)->default_value(0))
    ("calculateInterference", po::value<bool>(&calculateInterference)->default_value(0), "Calculate interference terms")
    ("fitFractionError", po::value<int>(&fitFractionError)->default_value(0),
     "Fit fraction errors are calculated by Monte-Carlo approach. 0 = disabled, >0 = MC precision");
    
    bool enablePlotting;
    bool plotDataSample;
    int plotSize;
    bool plotCorrectEfficiency;
    int plotNBins;
    bool plotAmplitude;
    bool gof_enable;
    int gof_mcPrecision;
    int gof_size;
    po::options_description config_plot("Plot settings");
    config_plot.add_options()
    ("plotting", po::value<bool>(&enablePlotting)->default_value(1), "Enable/Disable plotting")
    ("plotData", po::value<bool>(&plotDataSample)->default_value(1), "Enable/Disable plotting of data")
    ("plotSize", po::value<int>(&plotSize)->default_value(1000000), "Size of sample for amplitude plot")
    ("plotNBins", po::value<int>(&plotNBins)->default_value(100), "Number of bins")
    ("plotCorrectEfficiency", po::value<bool>(&plotCorrectEfficiency)->default_value(0), "Number of bins")
    ("plotAmplitude", po::value<bool>(&plotAmplitude)->default_value(1), "Enable/Disable plotting of amplitude")
    ("enableGoF", po::value<bool>(&gof_enable)->default_value(0), "Enable/Disable calculation of goodness-of-fit value")
    ("GoF_mcPrecision", po::value<int>(&gof_mcPrecision)->default_value(-1), "Number of MC events for PointToPoint GoF test calculation")
    ("GoF_size", po::value<int>(&gof_size)->default_value(10000), "Number of MC events for PointToPoint GoF test calculation");
    
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(config_inout);
    cmdline_options.add(config_eff).add(config_fit).add(config_plot);
    po::options_description config_file_options;
    config_file_options.add(config).add(config_inout).add(config_eff);
    config_file_options.add(config_fit).add(config_plot);
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
    
    po::positional_options_description p;
    store(po::command_line_parser(argc, argv).
          options(cmdline_options).positional(p).run(), vm);
    notify(vm);
    if (vm.count("help")) { cout << cmdline_options << "\n"; return 1; }
    
    ifstream ifs(config_file.c_str());
    if(ifs) {
        store(parse_config_file(ifs, config_file_options), vm);
        notify(vm);
    }
    std::string fileNamePrefix = outputDir+std::string("/")+outputFileName;
    std::string logFileName = fileNamePrefix+std::string(".log");
    if(disableFileLog) logFileName = "";
    Logging log(logFileName,boost::log::trivial::info); //initialize logging
    log.setLogLevel(logLevel);
    
    //check configuration
    assert(!outputDir.empty());
    if(trueModelFile.empty()){
        LOG(error) << "True model file not set";
    }
    if(fitModelFile.empty() && fittingMethod!="plotOnly"){
        LOG(error) << "Fit model file not set";
    }
    if(fittingMethod!="tree" && fittingMethod!="amplitude" && fittingMethod!="plotOnly") {
        LOG(error) << "Unknown fitting method: "
        <<fittingMethod; return 1;
    }
    
    LOG(info) << "Current path: " << getenv("PWD");
    po::notify(vm);

    
    if( ampMcPrecision == 0 ) ampMcPrecision = mcPrecision;
    
    PhysConst::createInstance();
    //initialize kinematics of decay
    Physics::DPKinematics::DalitzKinematics::createInstance("D0","K_S0","K-","K+");
    
    //initialize random generator
    std::shared_ptr<Generator> gen =
    std::shared_ptr<Generator>(new Physics::DPKinematics::RootGenerator(seed));
    
    RunManager run;
    run.setGenerator(gen);
    //======================= EFFICIENCY =============================
    auto eff = std::shared_ptr<Efficiency>(new UnitEfficiency());
    
    //Binned efficiency
    if(!efficiencyFile.empty()){
        TFile* tf= new TFile(TString(efficiencyFile));
        Physics::DPKinematics::SimpleEfficiency* effDalitz =
        (Physics::DPKinematics::SimpleEfficiency*) tf->Get(TString(efficiencyObject));
        tf->Close();
        eff = std::shared_ptr<Efficiency>(new Physics::DPKinematics::SimpleAngleEfficiency(effDalitz));
    }
    
    //Unbinned efficiency
    std::shared_ptr<Data> phspData, phspTrueData;
    if(!phspEfficiencyFile.empty()){
        //sample with accepted phsp events
        phspData =
        std::shared_ptr<Data>(
                              new DataReader::RootReader::RootReader(
                                                                     phspEfficiencyFile,
                                                                     phspEfficiencyFileTreeName,
                                                                     mcPrecision )
                              );
        phspData->reduceToPhsp();
    }
    //======================= AMPLITUDE =============================
    std::vector<std::shared_ptr<Data> > dataVec, trueDataVec;
    std::vector<std::shared_ptr<Amplitude> > ampVec, trueAmpVec;
    std::vector<double> fraction, trueFraction;
    
    // ========= SIGNAL AMPLITUDE ========
    //TRUE MODEL
    if(!trueModelFile.empty()){
        if(trueAmpOption.empty()) trueAmpOption = ampOption;
        createAmp("trueDzeroAmp",trueAmpVec,trueModelFile,eff,ampMcPrecision,trueAmpOption);
        if(trueAmpVec.size()==1) {
            trueDataVec.push_back( std::shared_ptr<Data>(new DataReader::RootReader::RootReader()) );
            trueFraction.push_back(trueSignalFraction);
        } else if(trueAmpVec.size()==2) {
            trueDataVec.push_back( std::shared_ptr<Data>(new DataReader::RootReader::RootReader()) );
            trueDataVec.push_back( std::shared_ptr<Data>(new DataReader::RootReader::RootReader()) );
            trueFraction.push_back(0.5*trueSignalFraction);
            trueFraction.push_back(0.5*trueSignalFraction);
        } else {
        }
    }
    
    //FIT MODEL
    if( !fitModelFile.empty() ) {
        createAmp("DzeroAmp",ampVec,fitModelFile,eff,ampMcPrecision,ampOption);
        if(ampVec.size()==1) {
            dataVec.push_back( std::shared_ptr<Data>(new DataReader::RootReader::RootReader()) );
            fraction.push_back(fitSignalFraction);
        } else if(ampVec.size()==2) {
            dataVec.push_back( std::shared_ptr<Data>(new DataReader::RootReader::RootReader()) );
            dataVec.push_back( std::shared_ptr<Data>(new DataReader::RootReader::RootReader()) );
            fraction.push_back(0.5*fitSignalFraction);
            fraction.push_back(0.5*fitSignalFraction);
        } else {
        }
    }
    
    // ========= BACKGROUND AMPLITUDE ========
    std::shared_ptr<Amplitude> trueBkgAmp, fitBkgAmp;
    //TRUE MODEL
    if(trueSignalFraction!=1. && !trueBkgModelFile.empty()){
        createAmp("trueBkg",trueAmpVec,trueBkgModelFile,eff,ampMcPrecision,"none");
        trueDataVec.push_back( std::shared_ptr<Data>(new DataReader::RootReader::RootReader()) );
        trueFraction.push_back(1-trueSignalFraction);
    } else if(trueSignalFraction!=1. && trueBkgModelFile.empty() && bkgFile.empty()){
        throw std::runtime_error("A true signal fraction was specified "
                                 "but no background model");
    } else {
        //No model specified. Use background sample
    }
    
    //FIT MODEL
    if(fitSignalFraction!=1.){
        if(fitBkgFile=="")
            throw std::runtime_error("A fit signal fraction was specified "
                                     "but no background model");
        
        createAmp("Bkg",ampVec,fitBkgFile,eff,ampMcPrecision,"none");
        dataVec.push_back( std::shared_ptr<Data>(new DataReader::RootReader::RootReader()) );
        fraction.push_back( 1-fitSignalFraction );
    }
    
    //======================= READING DATA =============================
    //sample is used for minimization
    std::shared_ptr<Data> sample( new DataReader::RootReader::RootReader() );
    
    //temporary samples for signal and background
    std::shared_ptr<Data> inputData, inputBkg;
    
    //smear total number of events
    if( smearNEvents && numEvents > 0 )
        numEvents += (int) gen->getGaussDist(0,sqrt(numEvents));
    
    if(trueSignalFraction!=1. && !bkgFile.empty()) {
        int numBkgEvents = (int)( (1-trueSignalFraction)*numEvents);
        LOG(info)<<"Reading background file...";
        std::shared_ptr<Data> inBkg(new DataReader::RootReader::RootReader(bkgFile,bkgFileTreeName));
        inBkg->reduceToPhsp();
        if(resetWeights) inBkg->resetWeights(); //resetting weights of requested
        inputBkg = inBkg->rndSubSet(numBkgEvents,gen);
        run.setBackground(inputBkg);
        inBkg = std::shared_ptr<Data>();
        sample->Add(*inputBkg);
    }
    if(!dataFile.empty()) {
        int numSignalEvents = numEvents;
        if( inputBkg ) numSignalEvents -= inputBkg->getNEvents();
        LOG(info)<<"Reading data file...";
        std::shared_ptr<Data> inD(new DataReader::RootReader::RootReader(dataFile,dataFileTreeName));
        inD->reduceToPhsp();
        if(resetWeights) inD->resetWeights(); //resetting weights of requested
        inputData = inD->rndSubSet(numSignalEvents,gen);
        inputData->setEfficiency(eff);
        run.setData(inputData);
        inD = std::shared_ptr<Data>();
        sample->Add(*inputData);
    }
    
    //======================== GENERATION =====================
    std::shared_ptr<Data> toyPhspData(new DataReader::RootReader::RootReader());//Toy sample
    run.setPhspSample(toyPhspData);
    run.generatePhsp(mcPrecision);
    toyPhspData->setEfficiency(eff);//set efficiency values for each event
    
    //generation with unbinned efficiency correction - use full sample
    if( !phspEfficiencyFile.empty() && (!inputData || ! inputBkg) ){
        //assume that efficiency of hit&miss is large 0.5%
        double phspSampleSize = numEvents/0.005;
        
        std::shared_ptr<Data> fullPhsp(
                                       new DataReader::RootReader::RootReader(
                                                                                      phspEfficiencyFile,
                                                                                      phspEfficiencyFileTreeName,
                                                                                      phspSampleSize
                                                                                      )
                                       );
        std::shared_ptr<Data> fullTruePhsp;
        if(!phspEfficiencyFileTrueTreeName.empty()){
            fullTruePhsp = std::shared_ptr<Data>(
                                                 new DataReader::RootReader::RootReader(
                                                                                                phspEfficiencyFile,
                                                                                                phspEfficiencyFileTrueTreeName,
                                                                                                phspSampleSize
                                                                                                )
                                                 );
            //Data* fullPhspRed = fullPhsp->EmptyClone();
            //Data* fullPhspTrueRed = fullTruePhsp->EmptyClone();
            //rndReduceSet(phspSampleSize,gen,fullPhsp.get(),fullPhspRed,
            //fullTruePhsp.get(),fullPhspTrueRed);
            //fullPhsp = std::shared_ptr<Data>(fullPhspRed);
            //fullTruePhsp = std::shared_ptr<Data>(fullPhspTrueRed);
        } else {
            //Data* fullPhspRed = fullPhsp->EmptyClone();
            //rndReduceSet(phspSampleSize,gen,fullPhsp.get(),fullPhspRed);
            //fullPhsp = std::shared_ptr<Data>(fullPhspRed);
        }
        run.setPhspSample( fullPhsp, fullTruePhsp );
    }
    if( !inputData ){
        LOG(info) <<"Generating samples for each amplitude!";
        run.SetAmplitudesData( trueAmpVec, trueFraction, trueDataVec );
        run.GenAmplitudesData( numEvents );
        for(int i=0; i<trueDataVec.size(); ++i){
            sample->Add(*trueDataVec.at(i));
        }
        LOG(info) <<"Sample size: "<<sample->getNEvents();
    }

    run.setPhspSample(std::shared_ptr<Data>()); //reset phsp sample to save memory
    
    sample->reduceToPhsp();
    
    LOG(info)<<"================== SETTINGS =================== ";
    
    LOG(info)<<"==== GENERAL";
    LOG(info)<<"Number of events: "<<numEvents;
    LOG(info)<<"Initial seed: "<<seed ;
    LOG(info)<<"Log level: "<<logLevel;
    LOG(info)<<"Output file: "<<fileNamePrefix
    << "-* [*.root,*.pdf,*.tex]";
    LOG(info)<<"==== INPUT/GENERATION";
    if(!dataFile.empty()){ //read in data
        LOG(info)<<"Data file: "<<dataFile
        <<" [tree="<<dataFileTreeName<<" ]";
    } else { //generate signal
        LOG(info)<<"true signal fraction: "<<trueSignalFraction;
        LOG(info)<<"True model file: "<<trueModelFile;
        LOG(info)<<"True amplitude option: "<<trueAmpOption;
    }
    if( trueSignalFraction != 1.0 ){
        if(!bkgFile.empty()){ //read in background
            LOG(info)<<"Background file: "<<bkgFile
            <<" [tree="<<bkgFileTreeName<<" ]";
        } else { //generate background
            LOG(info)<<"True background model file: "
            <<trueBkgModelFile;
        }
    }
    LOG(info)<<"Total events in input sample: "
    <<sample->getNEvents();
    
    LOG(info)<<"==== EFFICIENCY";
    if(!efficiencyFile.empty())
        LOG(info)<<"Binned Efficiency file: "
        <<efficiencyFile<<" [obj="<<efficiencyObject<<"]";
    if(!phspEfficiencyFile.empty()){
        LOG(info)<<"PHSP data input file: "<<phspEfficiencyFile
        <<" [tree="<<phspEfficiencyFileTreeName<<"]";
        if(!phspEfficiencyFileTrueTreeName.empty())
            LOG(info)<<"True tree name: "
            <<phspEfficiencyFileTrueTreeName;
        if(applySysCorrection)
            LOG(info)<<"Unbinned efficiency is corrected "
            "for tracking systematics!";
        LOG(info)<<"Unbinned efficiency correction is used!";
    }
    if(phspEfficiencyFile.empty() && !efficiencyFile.empty())
        LOG(info)<<"Binned efficiency correction is used!";
    if(!phspEfficiencyFile.empty() && !efficiencyFile.empty())
        LOG(info)<<"No efficiency correction!";
    
    LOG(info)<<"==== FIT";
    LOG(info)<<"MC precision: "<<mcPrecision;
    LOG(info)<<"Amplitude MC precision: "<<ampMcPrecision;
    LOG(info)<<"Fitting method: "<<fittingMethod;
    LOG(info)<<"Fit signal fraction: "<<fitSignalFraction;
    LOG(info)<<"Amplitude option: "<<ampOption;
    LOG(info)<<"Fit model file: "<<fitModelFile;
    if(!fitBkgFile.empty())
        LOG(info)<<"Fit bkg model file: "<<fitBkgFile;
    LOG(info)<<"Using Geneva as pre fitter: "<<usePreFitter;
    LOG(info)<<"Use MINOS: "<<useMinos;
    LOG(info)<<"Accurate errors on fit fractions: "
    <<fitFractionError;
    LOG(info)<<"Penalty scale: "<<penaltyScale;
    LOG(info)<<"==== GOF";
    LOG(info)<<"enable gof: "<<gof_enable;
    if(gof_enable){
        LOG(info)<<"size of sample: "<<gof_size;
        LOG(info)<<"mc precision: "<<gof_mcPrecision;
    }
    
    LOG(info)<<"==== PLOTTING";
    LOG(info)<<"enable plotting: "<<enablePlotting;
    if(enablePlotting){
        LOG(info)<<"plotting size: "<<plotSize;
        LOG(info)<<"number of bins: "<<plotNBins;
        LOG(info)<<"Correct samples for efficiency: "
        <<plotCorrectEfficiency;
        LOG(info)<<"plotDataSample: "<<plotDataSample;
        LOG(info)<<"plotAmplitude: "<<plotAmplitude;
    }
    LOG(info)<<"===============================================";
    
    std::shared_ptr<FitResult> result;
    ParameterList finalParList;
    std::shared_ptr<ControlParameter> esti;
    if( fittingMethod != "plotOnly") {
        //========================FITTING =====================
        ParameterList truePar, fitPar;
        Amplitude::FillAmpParameterToList(trueAmpVec,truePar);
        Amplitude::FillAmpParameterToList(ampVec,fitPar);
        
        //=== Constructing likelihood
        esti = std::shared_ptr<ControlParameter>(Estimator::MinLogLH::MinLogLH::createInstance(ampVec, fraction, sample, toyPhspData, phspData, 0, 0 ) );
        
        Estimator::MinLogLH::MinLogLH* contrPar = dynamic_cast<Estimator::MinLogLH::MinLogLH*>(&*(esti->Instance()));
        
        contrPar->setAmplitude(ampVec, fraction, sample, toyPhspData, phspData,
                               0, 0, false); //setting new trees in ControlParameter
        contrPar->setPenaltyScale(penaltyScale,0);
        
        if(fittingMethod == "tree") {
            contrPar->setUseFunctionTree(1);
            LOG(debug) << contrPar->getTree()->head()->to_str(20);
            //			auto trNode = contrPar->getTree()->head()->getChildNode("N");
            //			LOG(info) <<"Compare normalization: (asmndj)"
            //				<<" Tree: "<<trNode->getChildSingleValue("normFactor")
            //				<<" Amplitude(VEGAS): "<< ampVec.at(0)->GetNormalization();
        }
        
        if(useRandomStartValues) randomStartValues(fitPar);
        LOG(debug) << "Initial LH="
        << esti->controlParameter(fitPar) <<".";
        
        std::shared_ptr<Optimizer::Optimizer> preOpti, opti;
        std::shared_ptr<FitResult> preResult;
        if( usePreFitter ){
            LOG(info) << "Running Geneva as pre fitter!";
            setErrorOnParameterList(fitPar, 0.5, useMinos);
            //			preOpti = std::shared_ptr<Optimizer>(
            //					new GenevaIF(esti,dkskkDir+"/PWA/geneva-config/")
            //			);
            preResult = preOpti->exec(fitPar);
            preResult->setTrueParameters(truePar);
            preResult->print();
        }
        
        //Set start error of 0.05 for parameters
        setErrorOnParameterList(fitPar, 0.05, useMinos);
        
        auto minuitif = new Optimizer::Minuit2::MinuitIF(esti, fitPar);
        minuitif->SetHesse(useHesse);
        opti = std::shared_ptr<Optimizer::Optimizer>(minuitif);
        run.setOptimizer(opti);
        
        //====== STARTING MINIMIZATION ======
        result = run.startFit(fitPar);
        
        //====== FIT RESULT =======
        Optimizer::Minuit2::MinuitResult* minuitResult = dynamic_cast<Optimizer::Minuit2::MinuitResult*>(&*result);
        finalParList = result->getFinalParameters();
        Amplitude::UpdateAmpParameterList(esti->getAmplitudes(), finalParList);
        result->setTrueParameters(truePar);
        minuitResult->setUseCorrelatedErrors(fitFractionError);
        minuitResult->SetCalcInterference(calculateInterference);
        result->print();
        
        std::ofstream ofs(fileNamePrefix+std::string("-fitResult.xml"));
        boost::archive::xml_oarchive oa(ofs);
        oa << BOOST_SERIALIZATION_NVP(result);
        //ofs.close(); don't close! Otherwise file is not correctly written
        if(!disableFileLog){ // Write fit result
            //result->writeTeX(fileNamePrefix+std::string("-fitResult.tex"));
            
            //Save final amplitude
            Physics::AmplitudeSum::AmpSumIntensity* fitAmpSum =
            dynamic_cast<Physics::AmplitudeSum::AmpSumIntensity*>(&*ampVec.at(0));
            LOG(info) << "Average resonance width of fit model: "
            <<fitAmpSum->averageWidth();
            fitAmpSum->Save(fileNamePrefix+std::string("-Model.xml"));
        }
        
        if( !trueModelFile.empty() ){
            double confLevel = std::sqrt(
                                         2*( minuitResult->GetTrueLH() - minuitResult->GetFinalLH() )
                                         );
            LOG(info) << "CHI2 TEST: "
            <<" finalLH = "<<minuitResult->GetFinalLH()
            <<" trueLH = "<<minuitResult->GetTrueLH()
            <<" -> confLevel [sigma] = "<<confLevel;
        }
    } else if(!inputResult.empty()){
        LOG(info) << "Reading MinuitResult from "<<inputResult;
        std::shared_ptr<Optimizer::Minuit2::MinuitResult> inResult;
        std::ifstream ifs(inputResult);
        if(!ifs.good()) throw std::runtime_error("input stream not good!");
        boost::archive::xml_iarchive ia(ifs);
        ia >> BOOST_SERIALIZATION_NVP(inResult);
        ifs.close();
        
        inResult->SetAmplitude(ampVec);
        inResult->setUseCorrelatedErrors(fitFractionError);
        if(calculateInterference)
            inResult->SetCalcInterference(0);
        inResult->setFractions(ParameterList()); //reset fractions list
        result = inResult;
        result->print();
    }
    
    //Fill final parameters if minimization was not run
    if(!finalParList.GetNParameter()){
        ParameterList tmpList;
        Amplitude::FillAmpParameterToList(ampVec,tmpList);
        finalParList.DeepCopy(tmpList);
    }
    
    //======================= PLOTTING =============================
    if(enablePlotting){
        if(fittingMethod != "plotOnly")
            Amplitude::UpdateAmpParameterList(ampVec,finalParList);
        
        //------- phase-space sample
        std::shared_ptr<Data> pl_phspSample( new DataReader::RootReader::RootReader() );
        LOG(info) << "Plotting results...";
        if(!phspEfficiencyFile.empty()){ //unbinned plotting
            //sample with accepted phsp events
            pl_phspSample = std::shared_ptr<Data>(
                                                  new DataReader::RootReader::RootReader(
                                                                                                 phspEfficiencyFile,
                                                                                                 phspEfficiencyFileTreeName,
                                                                                                 plotSize
                                                                                                 )
                                                  );
            
            std::shared_ptr<Data> plotTruePhsp;
            if(!phspEfficiencyFileTrueTreeName.empty())
                plotTruePhsp = std::shared_ptr<Data>(
                                                     new DataReader::RootReader::RootReader(
                                                                                                    phspEfficiencyFile,
                                                                                                    phspEfficiencyFileTrueTreeName,
                                                                                                    plotSize
                                                                                                    )
                                                     );
            run.setPhspSample( pl_phspSample, plotTruePhsp );
            //make sure no efficiency is set
            Amplitude::SetAmpEfficiency(
                                        ampVec,
                                        std::shared_ptr<Efficiency>(new UnitEfficiency)
                                        );
        }else{ //binned plotting
            run.setPhspSample(pl_phspSample);
            run.generatePhsp(plotSize);//we generate a very large sample for plotting
        }
        //reduce sample to phsp
        pl_phspSample->reduceToPhsp();
        pl_phspSample->setEfficiency(eff);
        
        //-------- Instance of plotData
        plotData pl(fileNamePrefix,plotNBins);
        //set data sample
        pl.setData(sample);
        //set amplitude
        pl.setFitAmp(ampVec, fraction);
        //select components to plot
        if(fitBkgFile != "")
            pl.DrawComponent(ampVec.size()-1,kRed);
        pl.DrawComponent(ampVec.size()-2,kOrange);
        //		pl.setHitMissData(pl_sample);
        pl.setPhspData(pl_phspSample);
        
        pl.setCorrectEfficiency(plotCorrectEfficiency);
        //Fill histograms and create canvases
        pl.Plot();
    }
    LOG(info) << "FINISHED!";
    
    //Exit code is exit code of fit routine. 0 is good/ 1 is bad
    if(result) return result->hasFailed();
    return 0;
}

