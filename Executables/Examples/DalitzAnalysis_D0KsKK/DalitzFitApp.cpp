//! Fit application for D0->K_S0 K_ K- analysis
/*!
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

#include "boost/program_options.hpp"
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
using namespace ComPWA;
BOOST_CLASS_EXPORT(ComPWA::Optimizer::Minuit2::MinuitResult)
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Physics/DPKinematics/RootEfficiency.cpp"
#include "Physics/DPKinematics/SimpleEfficiency.hpp"
#include "Physics/DPKinematics/SimpleAngleEfficiency.hpp"
#include "Physics/DPKinematics/RootGenerator.hpp"

#include "PlotData.hpp"

using namespace std;
using namespace ComPWA;

//Expands '~' in file path to user directory
std::string expand_user(std::string p)
{
    std::string path = p;
    if (not path.empty() and path[0] == '~') {
        assert(path.size() == 1 or path[1] == '/');  // or other error handling
        char const* home = getenv("HOME");
        if( home || home==getenv("USERPROFILE") ) {
            path.replace(0, 1, home);
        } else {
            char const *hdrive = getenv("HOMEDRIVE"),
            *hpath = getenv("HOMEPATH");
            assert(hdrive);  // or other error handling
            assert(hpath);
            path.replace(0, 1, std::string(hdrive) + hpath);
        }
    }
    return path;
}

void createAmp(std::string name, std::vector<std::shared_ptr<Amplitude> >& ampV,
               std::string xmlInput, std::shared_ptr<Efficiency> eff,
               double mcPrecision, std::string ampOption)
{
    auto DzeroAmp = new ComPWA::Physics::AmplitudeSum::AmpSumIntensity(name,normStyle::one, eff, mcPrecision);
    DzeroAmp->Configure(expand_user(xmlInput));
    auto tmpAmp = std::shared_ptr<ComPWA::Physics::AmplitudeSum::AmpSumIntensity>(DzeroAmp);
    ampV.push_back(tmpAmp);
}

void setErrorOnParameterList(ParameterList& list, double error, bool asym)
{
    for(unsigned int i=0; i<list.GetNDouble(); i++){
        std::shared_ptr<DoubleParameter> p = list.GetDoubleParameter(i);
        if(p->IsFixed()) {
            p->SetError(0.0);
            continue;
        }
        if(asym)
            list.GetDoubleParameter(i)->SetError(
                                                 error,error
                                                 );
        else list.GetDoubleParameter(i)->SetError(error);
    }
}

void randomStartValues(ParameterList& fitPar)
{
    for(unsigned int i=0; i<fitPar.GetNDouble(); i++){
        std::shared_ptr<DoubleParameter> p = fitPar.GetDoubleParameter(i);
        if(p->IsFixed()) continue;
        double min = -999, max = 999;
        if(p->UseBounds()){
            min = p->GetMinValue();
            max = p->GetMaxValue();
        }
        p->SetValue(gRandom->Uniform(min,max));
    }
    std::cout<<"Randomizing parameter list. New list:"<<fitPar<<std::endl;
    return;
}

/*****************************************************************************/
/**
 * The main function.
 */
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
	("nEvents", po::value<int>(&numEvents)->default_value(-1), "set of events per fit")
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
	("trueModelFile", po::value<std::string>(&trueModelFile)->default_value(""), "set XML model file of true distribution")
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
	("fitModelFile", po::value<std::string>(&fitModelFile)->default_value(""), "Set XML model file of fit distribution")
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

	std::string dkskkDir = getenv("DKSKK_DIR");


	//check configuration
	assert(!outputDir.empty());
	if(trueModelFile.empty()){
		BOOST_LOG_TRIVIAL(error) << "True model file not set";
	}
	if(fitModelFile.empty() && fittingMethod!="plotOnly"){
		BOOST_LOG_TRIVIAL(error) << "Fit model file not set";
	}
	if(fittingMethod!="tree" && fittingMethod!="amplitude" && fittingMethod!="plotOnly") {
		BOOST_LOG_TRIVIAL(error) << "Unknown fitting method: "
				<<fittingMethod; return 1;
	}
	if(logLevel!="info"&&logLevel!="warning"&&logLevel!="error"&&logLevel!="debug"){
		BOOST_LOG_TRIVIAL(error) << "Unknown log level: "<<logLevel; return 1;
	}
	boost::log::trivial::severity_level logLv;
	if(logLevel == "info") logLv = boost::log::trivial::info;
	if(logLevel == "error") logLv = boost::log::trivial::error;
	if(logLevel == "debug") logLv = boost::log::trivial::debug;
	if(logLevel == "warning") logLv = boost::log::trivial::warning;
	log.setLogLevel(logLv);

	po::notify(vm);

	//initialize kinematics of decay
	ComPWA::Physics::DPKinematics::DalitzKinematics::createInstance("D0","K_S0","K-","K+");

	if( ampMcPrecision == 0 ) ampMcPrecision = mcPrecision;

	//initialize random generator
	std::shared_ptr<Generator> gen =
			std::shared_ptr<Generator>(new ComPWA::Physics::DPKinematics::RootGenerator(seed));

	RunManager run;
	run.setGenerator(gen);

	//======================= EFFICIENCY =============================
	std::shared_ptr<Efficiency> eff =
			std::shared_ptr<Efficiency>(new UnitEfficiency());;

	//Binned efficiency
	if(!efficiencyFile.empty()){
		TFile* tf= new TFile(TString(efficiencyFile));
		ComPWA::Physics::DPKinematics::SimpleEfficiency* effDalitz =
				(ComPWA::Physics::DPKinematics::SimpleEfficiency*) tf->Get(TString(efficiencyObject));
		tf->Close();
		eff = std::shared_ptr<Efficiency>(new ComPWA::Physics::DPKinematics::SimpleAngleEfficiency(effDalitz));
	}

	//Unbinned efficiency
	std::shared_ptr<Data> phspData, phspTrueData;
	if(!phspEfficiencyFile.empty()){
		//sample with accepted phsp events
		phspData = std::shared_ptr<Data>(
				new ComPWA::DataReader::RootReader::RootReader(
						phspEfficiencyFile,
						phspEfficiencyFileTreeName,
						mcPrecision
				)
		);
		phspData->reduceToPhsp();
		/* Set total efficiency of phsp sample. Necessary for consitency check
		 * with unbinned correction. Doesn't influence the fit result. */
		//phspData->resetEfficiency( 1.0 );
		phspData->resetEfficiency( 0.0135554 );
//		if(applySysCorrection){
//			MomentumCorrection* trkSys = getTrackingCorrection();
//			MomentumCorrection* pidSys = getPidCorrection();
//			trkSys->Print();
//			pidSys->Print();
//			phspData->applyCorrection(*trkSys);
//			phspData->applyCorrection(*pidSys);
//		}
		if( fittingMethod == "amplitude" ){
			BOOST_LOG_TRIVIAL(error)<<"ATTENTION: amplitude fit with unbinned "
					"efficiency correction needs UnitEfficiency at this point! "
					"During plotting the normalization is not correct";
		}
	}
	//======================= AMPLITUDE =============================
	std::vector<std::shared_ptr<Data> > dataVec, trueDataVec;
	std::vector<std::shared_ptr<Amplitude> > ampVec, trueAmpVec;
	std::vector<double> fraction, trueFraction;

	double trueNorm = -1, fitNorm = -1, trueBkgNorm = -1, fitBkgNorm = -1;
	// ========= SIGNAL AMPLITUDE ========
	//TRUE MODEL
	if(!trueModelFile.empty()){
		if(trueAmpOption.empty()) trueAmpOption = ampOption;
		createAmp("trueDzeroAmp",trueAmpVec,trueModelFile,eff,ampMcPrecision,trueAmpOption);
		if(trueAmpVec.size()==1) {
			trueDataVec.push_back( std::shared_ptr<Data>(new ComPWA::DataReader::RootReader::RootReader()) );
			trueFraction.push_back(trueSignalFraction);
		} else if(trueAmpVec.size()==2) {
			trueDataVec.push_back( std::shared_ptr<Data>(new ComPWA::DataReader::RootReader::RootReader()) );
			trueDataVec.push_back( std::shared_ptr<Data>(new ComPWA::DataReader::RootReader::RootReader()) );
			trueFraction.push_back(0.5*trueSignalFraction);
			trueFraction.push_back(0.5*trueSignalFraction);
		} else {
		}
	}

	//FIT MODEL
	if( !fitModelFile.empty() ) {
		createAmp("DzeroAmp",ampVec,fitModelFile,eff,ampMcPrecision,ampOption);
		if(ampVec.size()==1) {
			dataVec.push_back( std::shared_ptr<Data>(new ComPWA::DataReader::RootReader::RootReader()) );
			fraction.push_back(fitSignalFraction);
		} else if(ampVec.size()==2) {
			dataVec.push_back( std::shared_ptr<Data>(new ComPWA::DataReader::RootReader::RootReader()) );
			dataVec.push_back( std::shared_ptr<Data>(new ComPWA::DataReader::RootReader::RootReader()) );
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
		trueDataVec.push_back( std::shared_ptr<Data>(new ComPWA::DataReader::RootReader::RootReader()) );
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
		dataVec.push_back( std::shared_ptr<Data>(new ComPWA::DataReader::RootReader::RootReader()) );
		fraction.push_back( 1-fitSignalFraction );
	}

	//======================= READING DATA =============================
	//sample is used for minimization
	std::shared_ptr<Data> sample( new ComPWA::DataReader::RootReader::RootReader() );

	//temporary samples for signal and background
	std::shared_ptr<Data> inputData, inputBkg;

	//smear total number of events
	if( smearNEvents && numEvents > 0 )
		numEvents += (int) gen->getGaussDist(0,sqrt(numEvents));

	if(trueSignalFraction!=1. && !bkgFile.empty()) {
		int numBkgEvents = (int)( (1-trueSignalFraction)*numEvents);
		BOOST_LOG_TRIVIAL(info)<<"Reading background file...";
		std::shared_ptr<Data> inBkg(new ComPWA::DataReader::RootReader::RootReader(bkgFile,bkgFileTreeName));
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
		BOOST_LOG_TRIVIAL(info)<<"Reading data file...";
		std::shared_ptr<Data> inD(new ComPWA::DataReader::RootReader::RootReader(dataFile,dataFileTreeName));
		inD->reduceToPhsp();
		if(resetWeights) inD->resetWeights(); //resetting weights of requested
		inputData = inD->rndSubSet(numSignalEvents,gen);
		inputData->setEfficiency(eff);
		run.setData(inputData);
		inD = std::shared_ptr<Data>();
		sample->Add(*inputData);
	}

	//======================== GENERATION =====================
	std::shared_ptr<Data> toyPhspData(new ComPWA::DataReader::RootReader::RootReader());//Toy sample
	run.setPhspSample(toyPhspData);
	run.generatePhsp(mcPrecision);
	toyPhspData->setEfficiency(eff);//set efficiency values for each event

	//generation with unbinned efficiency correction - use full sample
	if( !phspEfficiencyFile.empty() && (!inputData || ! inputBkg) ){
		//assume that efficiency of hit&miss is large 0.5%
		double phspSampleSize = numEvents/0.005;

		std::shared_ptr<Data> fullPhsp(
				new ComPWA::DataReader::RootReader::RootReader(
					phspEfficiencyFile,
					phspEfficiencyFileTreeName,
					phspSampleSize
					)
				);
		std::shared_ptr<Data> fullTruePhsp;
		if(!phspEfficiencyFileTrueTreeName.empty()){
			fullTruePhsp = std::shared_ptr<Data>(
					new ComPWA::DataReader::RootReader::RootReader(
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
		BOOST_LOG_TRIVIAL(info) <<"Generating samples for each amplitude!";
		run.SetAmplitudesData( trueAmpVec, trueFraction, trueDataVec );
		run.GenAmplitudesData( numEvents );
		for(int i=0; i<trueDataVec.size(); ++i){
			sample->Add(*trueDataVec.at(i));
		}
		BOOST_LOG_TRIVIAL(info) <<"Sample size: "<<sample->getNEvents();
	}
	//generate signal events
	//	if( !inputData && numSignalEvents ) {//generate signal from model
	//		if(trueModelFile.empty())
	//			throw std::runtime_error(
	//					"No true signal model given! Can't generate toy MC!");
	//		inputData = std::shared_ptr<Data>( new RootReader() );
	//		run.setData(inputData);
	//
	//		/*WORKAROUND:
	//		 * If toyPhspData is too small (mcPrecision too small),
	//		 * comment this out to generate on the fly */
	//		//run.setPhspSample(std::shared_ptr<Data>());
	//		run.generate( numSignalEvents  );
	//		//run.setPhspSample(toyPhspData);
	//	}
	//	//generate background events
	//	if( !inputBkg && numBkgEvents ){//generate background from model
	//		if(trueBkgFile.empty())
	//			throw std::runtime_error(
	//					"No true background model given! Can't generate background toy MC!");
	//		inputBkg= std::shared_ptr<Data>( new RootReader() );
	//		run.setBackground(inputBkg);
	//		run.generateBkg(numBkgEvents);
	//	}
	run.setPhspSample(std::shared_ptr<Data>()); //reset phsp sample to save memory

	//add signal and background sample
	sample->reduceToPhsp();

	BOOST_LOG_TRIVIAL(info)<<"================== SETTINGS =================== ";

	BOOST_LOG_TRIVIAL(info)<<"==== GENERAL";
	BOOST_LOG_TRIVIAL(info)<<"Number of events: "<<numEvents;
	BOOST_LOG_TRIVIAL(info)<<"Initial seed: "<<seed ;
	BOOST_LOG_TRIVIAL(info)<<"Log level: "<<logLevel;
	BOOST_LOG_TRIVIAL(info)<<"Output file: "<<fileNamePrefix
		<< "-* [*.root,*.pdf,*.tex]";
	BOOST_LOG_TRIVIAL(info)<<"==== INPUT/GENERATION";
	if(!dataFile.empty()){ //read in data
		BOOST_LOG_TRIVIAL(info)<<"Data file: "<<dataFile
			<<" [tree="<<dataFileTreeName<<" ]";
	} else { //generate signal
		BOOST_LOG_TRIVIAL(info)<<"true signal fraction: "<<trueSignalFraction;
		BOOST_LOG_TRIVIAL(info)<<"True model file: "<<trueModelFile;
		BOOST_LOG_TRIVIAL(info)<<"True amplitude option: "<<trueAmpOption;
	}
	if( trueSignalFraction != 1.0 ){
		if(!bkgFile.empty()){ //read in background
			BOOST_LOG_TRIVIAL(info)<<"Background file: "<<bkgFile
				<<" [tree="<<bkgFileTreeName<<" ]";
		} else { //generate background
			BOOST_LOG_TRIVIAL(info)<<"True background model file: "
				<<trueBkgModelFile;
		}
	}
	BOOST_LOG_TRIVIAL(info)<<"Total events in input sample: "
		<<sample->getNEvents();

	BOOST_LOG_TRIVIAL(info)<<"==== EFFICIENCY";
	if(!efficiencyFile.empty())
		BOOST_LOG_TRIVIAL(info)<<"Binned Efficiency file: "
			<<efficiencyFile<<" [obj="<<efficiencyObject<<"]";
	if(!phspEfficiencyFile.empty()){
		BOOST_LOG_TRIVIAL(info)<<"PHSP data input file: "<<phspEfficiencyFile
			<<" [tree="<<phspEfficiencyFileTreeName<<"]";
		if(!phspEfficiencyFileTrueTreeName.empty())
			BOOST_LOG_TRIVIAL(info)<<"True tree name: "
				<<phspEfficiencyFileTrueTreeName;
		if(applySysCorrection)
			BOOST_LOG_TRIVIAL(info)<<"Unbinned efficiency is corrected "
				"for tracking systematics!";
		BOOST_LOG_TRIVIAL(info)<<"Unbinned efficiency correction is used!";
	}
	if(phspEfficiencyFile.empty() && !efficiencyFile.empty())
		BOOST_LOG_TRIVIAL(info)<<"Binned efficiency correction is used!";
	if(!phspEfficiencyFile.empty() && !efficiencyFile.empty())
		BOOST_LOG_TRIVIAL(info)<<"No efficiency correction!";

	BOOST_LOG_TRIVIAL(info)<<"==== FIT";
	BOOST_LOG_TRIVIAL(info)<<"MC precision: "<<mcPrecision;
	BOOST_LOG_TRIVIAL(info)<<"Amplitude MC precision: "<<ampMcPrecision;
	BOOST_LOG_TRIVIAL(info)<<"Fitting method: "<<fittingMethod;
	BOOST_LOG_TRIVIAL(info)<<"Fit signal fraction: "<<fitSignalFraction;
	BOOST_LOG_TRIVIAL(info)<<"Amplitude option: "<<ampOption;
	BOOST_LOG_TRIVIAL(info)<<"Fit model file: "<<fitModelFile;
	if(!fitBkgFile.empty())
		BOOST_LOG_TRIVIAL(info)<<"Fit bkg model file: "<<fitBkgFile;
	BOOST_LOG_TRIVIAL(info)<<"Using Geneva as pre fitter: "<<usePreFitter;
	BOOST_LOG_TRIVIAL(info)<<"Use MINOS: "<<useMinos;
	BOOST_LOG_TRIVIAL(info)<<"Accurate errors on fit fractions: "
		<<fitFractionError;
	BOOST_LOG_TRIVIAL(info)<<"Penalty scale: "<<penaltyScale;
	BOOST_LOG_TRIVIAL(info)<<"==== GOF";
	BOOST_LOG_TRIVIAL(info)<<"enable gof: "<<gof_enable;
	if(gof_enable){
		BOOST_LOG_TRIVIAL(info)<<"size of sample: "<<gof_size;
		BOOST_LOG_TRIVIAL(info)<<"mc precision: "<<gof_mcPrecision;
	}

	BOOST_LOG_TRIVIAL(info)<<"==== PLOTTING";
	BOOST_LOG_TRIVIAL(info)<<"enable plotting: "<<enablePlotting;
	if(enablePlotting){
		BOOST_LOG_TRIVIAL(info)<<"plotting size: "<<plotSize;
		BOOST_LOG_TRIVIAL(info)<<"number of bins: "<<plotNBins;
		BOOST_LOG_TRIVIAL(info)<<"Correct samples for efficiency: "
				<<plotCorrectEfficiency;
		BOOST_LOG_TRIVIAL(info)<<"plotDataSample: "<<plotDataSample;
		BOOST_LOG_TRIVIAL(info)<<"plotAmplitude: "<<plotAmplitude;
	}
	BOOST_LOG_TRIVIAL(info)<<"===============================================";

	std::shared_ptr<FitResult> result;
	ParameterList finalParList;
	std::shared_ptr<ComPWA::ControlParameter> esti;
	if( fittingMethod != "plotOnly") {
		//========================FITTING =====================
		ParameterList truePar, fitPar;
		Amplitude::FillAmpParameterToList(trueAmpVec,truePar);
		Amplitude::FillAmpParameterToList(ampVec,fitPar);

		//=== Constructing likelihood
		esti = std::shared_ptr<ComPWA::ControlParameter>(                                                         ComPWA::Estimator::MinLogLH::MinLogLH::createInstance(
					ampVec, fraction, sample, toyPhspData, phspData, 0, 0
					)
				);

        ComPWA::Estimator::MinLogLH::MinLogLH* contrPar = dynamic_cast<ComPWA::Estimator::MinLogLH::MinLogLH*>(&*(esti->Instance()));

		contrPar->setAmplitude(ampVec, fraction, sample, toyPhspData, phspData,
				0, 0, false); //setting new trees in ControlParameter
		contrPar->setPenaltyScale(penaltyScale,0);

		if(fittingMethod == "tree") {
			contrPar->setUseFunctionTree(1);
			BOOST_LOG_TRIVIAL(debug) << contrPar->getTree()->head()->to_str(20);
//			auto trNode = contrPar->getTree()->head()->getChildNode("N");
//			BOOST_LOG_TRIVIAL(info) <<"Compare normalization: (asmndj)"
//				<<" Tree: "<<trNode->getChildSingleValue("normFactor")
//				<<" Amplitude(VEGAS): "<< ampVec.at(0)->GetNormalization();
		}

		if(useRandomStartValues) randomStartValues(fitPar);
		BOOST_LOG_TRIVIAL(debug) << "Initial LH="
			<< esti->controlParameter(fitPar) <<".";

        std::shared_ptr<ComPWA::Optimizer::Optimizer> preOpti, opti;
		std::shared_ptr<FitResult> preResult;
		if( usePreFitter ){
			BOOST_LOG_TRIVIAL(info) << "Running Geneva as pre fitter!";
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

        auto minuitif = new ComPWA::Optimizer::Minuit2::MinuitIF(esti, fitPar);
		minuitif->SetHesse(useHesse);
        opti = std::shared_ptr<ComPWA::Optimizer::Optimizer>(minuitif);
		run.setOptimizer(opti);

		//====== STARTING MINIMIZATION ======
		result = run.startFit(fitPar);

		//====== FIT RESULT =======
        ComPWA::Optimizer::Minuit2::MinuitResult* minuitResult = dynamic_cast<ComPWA::Optimizer::Minuit2::MinuitResult*>(&*result);
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
			ComPWA::Physics::AmplitudeSum::AmpSumIntensity* fitAmpSum =
				dynamic_cast<ComPWA::Physics::AmplitudeSum::AmpSumIntensity*>(&*ampVec.at(0));
			BOOST_LOG_TRIVIAL(info) << "Average resonance width of fit model: "
				<<fitAmpSum->averageWidth();
			fitAmpSum->Save(fileNamePrefix+std::string("-Model.xml"));
		}

		if( !trueModelFile.empty() ){
			double confLevel = std::sqrt(
					2*( minuitResult->GetTrueLH() - minuitResult->GetFinalLH() )
			);
			BOOST_LOG_TRIVIAL(info) << "CHI2 TEST: "
					<<" finalLH = "<<minuitResult->GetFinalLH()
					<<" trueLH = "<<minuitResult->GetTrueLH()
					<<" -> confLevel [sigma] = "<<confLevel;
		}
	} else if(!inputResult.empty()){
		BOOST_LOG_TRIVIAL(info) << "Reading MinuitResult from "<<inputResult;
		std::shared_ptr<ComPWA::Optimizer::Minuit2::MinuitResult> inResult;
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
		std::shared_ptr<Data> pl_phspSample( new ComPWA::DataReader::RootReader::RootReader() );
		BOOST_LOG_TRIVIAL(info) << "Plotting results...";
		if(!phspEfficiencyFile.empty()){ //unbinned plotting
			//sample with accepted phsp events
			pl_phspSample = std::shared_ptr<Data>(
					new ComPWA::DataReader::RootReader::RootReader(
						phspEfficiencyFile,
						phspEfficiencyFileTreeName,
						plotSize
						)
					);

			std::shared_ptr<Data> plotTruePhsp;
			if(!phspEfficiencyFileTrueTreeName.empty())
				plotTruePhsp = std::shared_ptr<Data>(
						new ComPWA::DataReader::RootReader::RootReader(
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

		//-------- Hit&miss sample for plotting
		//		std::vector<std::shared_ptr<Data> > pl_sampleVec;
		//		std::shared_ptr<Data> pl_sample( new RootReader() );
		//		for(int i=0; i<ampVec.size(); ++i)
		//			pl_sampleVec.push_back( std::shared_ptr<Data>( new RootReader() ) );
		//		run.SetAmplitudesData( ampVec, fraction, pl_sampleVec ) ;
		//		run.GenAmplitudesData( 12000 );
		//		//merge samples
		//		for(int i=0; i<pl_sampleVec.size(); ++i){
		//			pl_sample->Add(*pl_sampleVec.at(i));
		//		}
		//		//reduce sample to phsp
		//		pl_sample->reduceToPhsp();
		//		pl_sample->setEfficiency(eff);

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
	BOOST_LOG_TRIVIAL(info) << "FINISHED!";

	//Exit code is exit code of fit routine. 0 is good/ 1 is bad
	if(result) return result->hasFailed();
	return 0;
}

