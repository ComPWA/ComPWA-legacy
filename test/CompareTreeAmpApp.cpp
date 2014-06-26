//-------------------------------------------------------------------------------
// Copyright (c) 2014 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff
//-------------------------------------------------------------------------------
/** Test application to validate tree and amplitude
 * @file CompareTreeAmpAmpp.cpp
 * Application to validate that fitting using the function tree and using the amplitude produces the
 * same result. The decay D0->K_S0 K+ K- is used as example and the Breit-Wigner description of the
 * phi(1020) and the Flatte description of the a_0(980)0 are compared between tree and amplitude.
 */

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

#include "boost/program_options.hpp"

// Root header files go here
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH2D.h"

//Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/RunManager.hpp"
#include "Core/Efficiency.hpp"
#include "Core/FunctionTree.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Physics/DPKinematics/RootEfficiency.cpp"
#include "Physics/DPKinematics/RootGenerator.cpp"
#include "Core/TableFormater.hpp"
#include "Core/AbsParameter.hpp"
#include "Core/Logging.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
	int seed = 3041; //default seed

	unsigned int mcPrecision = 100000; //number of calls for numeric integration and number of events for phsp integration
	Logging log("log-compareTreeAmp.txt",boost::log::trivial::info); //initialize logging
	//initialize kinematics of decay
	DalitzKinematics::createInstance("D0","K_S0","K-","K+");//setup kinematics
	//initialize random generator
	std::shared_ptr<Generator> gen = std::shared_ptr<Generator>(new RootGenerator(seed));//initialize random generator

	RunManager run;
	run.setGenerator(gen);
	//======================= DATA =============================
	allMasses myToyPhspMasses, myEvtMasses;

	unsigned int numEvents = 300;//data size to be generated
	std::shared_ptr<Data> inputData(new RootReader("out.root", false,"data",false)); //empty file: run generation before fit

	//======================= EFFICIENCY =============================
	std::shared_ptr<Efficiency> eff(new UnitEfficiency());

	//======================= AMPLITUDE =============================
	//true amplitude model
	std::string trueModelFile = "test/CompareTreeAmp-model.xml";
	AmplitudeSetup iniTrue(trueModelFile);//put start parameters here
	std::shared_ptr<Amplitude> trueAmp( new AmpSumIntensity(iniTrue, AmpSumIntensity::normStyle::one, eff, mcPrecision) );
	//fit amplitude model
	std::string fitModelFile = trueModelFile;
	AmplitudeSetup ini(fitModelFile);//put start parameters here
	AmplitudeSetup iniTree(fitModelFile);//put start parameters here
	AmpSumIntensity* fitAmpPtr = new AmpSumIntensity(ini, AmpSumIntensity::normStyle::one, eff, mcPrecision);
	std::shared_ptr<Amplitude> fitAmp(fitAmpPtr);
	AmpSumIntensity* fitAmpTreePtr = new AmpSumIntensity(iniTree, AmpSumIntensity::normStyle::one, eff, mcPrecision);
	std::shared_ptr<Amplitude> fitAmpTree(fitAmpTreePtr);

	run.setAmplitude(trueAmp);//set true model here for generation
	std::shared_ptr<Data> toyPhspData(new RootReader("out.root", false,"mc",false));//empty phsp sample
	run.setPhspSample(toyPhspData);
	if( !toyPhspData->getNEvents() ) {
		run.generatePhsp(mcPrecision);
		myToyPhspMasses = allMasses(toyPhspData->getMasses());
		myToyPhspMasses.setEfficiency(eff); //set efficiency for all datapoints
	}

	run.setData(inputData);
	if( !inputData->getNEvents() ) {
		run.generate(numEvents);
		myEvtMasses = allMasses(inputData->getMasses());
	}

	//setup trees
	std::shared_ptr<FunctionTree> physicsTree = fitAmpTree->functionTree(myEvtMasses,myToyPhspMasses);
	if(!physicsTree){
		BOOST_LOG_TRIVIAL(error)<<"Physics Trees not setup correctly, quitting";
		return 0;
	}
	//	std::shared_ptr<FunctionTree> phspTree = fitAmpTree->phspTree(myPhspMasses,myToyPhspMasses);//unbinned efficiency
	std::shared_ptr<FunctionTree> phspTree = fitAmpTree->phspTree(myToyPhspMasses); //binned efficiency
	if(!phspTree){
		BOOST_LOG_TRIVIAL(error)<<"Phsp Trees not setup correctly, quitting";
		return 0;
	}
	BOOST_LOG_TRIVIAL(debug)<<"Check Trees: ";
	if(!physicsTree->sanityCheck()) return 0;
	if(!phspTree->sanityCheck()) return 0;

	//======================= PARAMETERS =============================
	ParameterList fitPar;
	fitAmp->fillStartParVec(fitPar);
	ParameterList fitParTree;
	fitAmpTree->fillStartParVec(fitParTree);
	ParameterList truePar;
	trueAmp->fillStartParVec(truePar); //true values

	for(unsigned int i=0; i<fitParTree.GetNDouble(); i++){
		fitParTree.GetDoubleParameter(i)->SetError( std::shared_ptr<ParError<double>>(new SymError<double>(.1)) );
	}

	fitParTree.GetDoubleParameter("mag_a_0(980)0")->FixParameter(1);
	fitParTree.GetDoubleParameter("phase_a_0(980)0")->FixParameter(1);

	fitParTree.GetDoubleParameter("mass_a_0(980)0")->FixParameter(1);
	fitParTree.GetDoubleParameter("width_a_0(980)0")->FixParameter(1);
	fitParTree.GetDoubleParameter("mass_a_0(980)+")->FixParameter(1);
	fitParTree.GetDoubleParameter("width_a_0(980)+")->FixParameter(1);
	fitParTree.GetDoubleParameter("mass_phi(1020)")->FixParameter(1);
	fitParTree.GetDoubleParameter("width_phi(1020)")->FixParameter(1);
	fitParTree.GetDoubleParameter("mass_f_0(1400)")->FixParameter(1);
	fitParTree.GetDoubleParameter("width_f_0(1400)")->FixParameter(1);

	for(unsigned int i=0; i<fitPar.GetNDouble(); i++){
		fitPar.GetDoubleParameter(i)->SetError( std::shared_ptr<ParError<double>>(new SymError<double>(.1)) );
	}

	fitPar.GetDoubleParameter("mag_a_0(980)0")->FixParameter(1);
	fitPar.GetDoubleParameter("phase_a_0(980)0")->FixParameter(1);

	fitPar.GetDoubleParameter("mass_a_0(980)0")->FixParameter(1);
	fitPar.GetDoubleParameter("width_a_0(980)0")->FixParameter(1);
	fitPar.GetDoubleParameter("mass_a_0(980)+")->FixParameter(1);
	fitPar.GetDoubleParameter("width_a_0(980)+")->FixParameter(1);
	fitPar.GetDoubleParameter("mass_phi(1020)")->FixParameter(1);
	fitPar.GetDoubleParameter("width_phi(1020)")->FixParameter(1);
	fitPar.GetDoubleParameter("mass_f_0(1400)")->FixParameter(1);
	fitPar.GetDoubleParameter("width_f_0(1400)")->FixParameter(1);

	fitAmpTree->setParameterList(fitParTree);
	fitAmp->setParameterList(fitPar);
//	std::cout<<fitPar<<std::endl;
	fitAmp->printAmps();

	BOOST_LOG_TRIVIAL(info)<<"Entries in data file: "<<inputData->getNEvents();
	BOOST_LOG_TRIVIAL(info)<<"True model file: "<<trueModelFile ;
	BOOST_LOG_TRIVIAL(info)<<"Fit model file: "<<fitModelFile ;

	//======================= TREE FIT =============================
	std::shared_ptr<ControlParameter> esti(MinLogLH::createInstance(physicsTree, phspTree, myEvtMasses.nEvents));
	double initialLHTree = esti->controlParameter(fitParTree);
	std::shared_ptr<Optimizer> optiTree(new MinuitIF(esti, fitParTree));
	run.setOptimizer(optiTree);

	BOOST_LOG_TRIVIAL(info)<<physicsTree<<std::endl;
	BOOST_LOG_TRIVIAL(info)<<phspTree<<std::endl;


	//======================= Compare tree and amplitude =============================
	std::shared_ptr<AmpRelBreitWignerRes> phiRes = std::dynamic_pointer_cast<AmpRelBreitWignerRes>(fitAmpPtr->getResonance("phi(1020)"));
	double phimag = fitAmpPtr->getMagnitude("phi(1020)");
	double phiphase = fitAmpPtr->getPhase("phi(1020)");
	std::complex<double> phiCoeff(phimag*cos(phiphase),phimag*sin(phiphase));
	std::shared_ptr<AmpFlatteRes> a0Res = std::dynamic_pointer_cast<AmpFlatteRes>(fitAmpPtr->getResonance("a_0(980)0"));
	double a0mag = fitAmpPtr->getMagnitude("a_0(980)0");
	double a0phase = fitAmpPtr->getPhase("a_0(980)0");
	std::complex<double> a0Coeff(a0mag*cos(a0phase),a0mag*sin(a0phase));

	dataPoint point(inputData->getEvent(0)); //first datapoint in sample
	ParameterList intens = fitAmpPtr->intensity(point);

	BOOST_LOG_TRIVIAL(info) <<"===========================================";
	BOOST_LOG_TRIVIAL(info) <<"Compare values: (use first event of data sample) TREE/AMPLITUDE";
	BOOST_LOG_TRIVIAL(info) <<"===========================================";
	BOOST_LOG_TRIVIAL(info) <<"Intensity: "<<physicsTree->head()->getChildValue("Intens").real()
			<<"/"<<*intens.GetDoubleParameter(0);
	BOOST_LOG_TRIVIAL(info) <<"================= phi(1020) ==========================";
	BOOST_LOG_TRIVIAL(info) <<"Reso_phi(1020): "<<physicsTree->head()->getChildValue("Reso_phi(1020)")
			<<"/"<<phiRes->evaluate(point)*phiCoeff;
	BOOST_LOG_TRIVIAL(info) <<"BW_phi(1020): "<<physicsTree->head()->getChildValue("BW_phi(1020)")
			<<"/"<<phiRes->evaluateAmp(point)*phiRes->GetNormalization();
	BOOST_LOG_TRIVIAL(info) <<"N_phi(1020): "<<physicsTree->head()->getChildValue("N_phi(1020)").real()
			<<"/"<<phiRes->GetNormalization();
	BOOST_LOG_TRIVIAL(info) <<"RelBW_phi(1020): "<<physicsTree->head()->getChildValue("RelBW_phi(1020)")
			<<"/"<<phiRes->evaluateAmp(point);
	BOOST_LOG_TRIVIAL(info) <<"AngD_phi(1020): "<<physicsTree->head()->getChildValue("AngD_phi(1020)").real()
			<<"/"<<phiRes->evaluateWignerD(point);
	BOOST_LOG_TRIVIAL(info) <<"================= a_0(980)0 ==========================";
	BOOST_LOG_TRIVIAL(info) <<"Reso_a_0(980)0: "<<physicsTree->head()->getChildValue("Reso_a_0(980)0")
			<<"/"<<a0Res->evaluate(point)*a0Coeff;
	BOOST_LOG_TRIVIAL(info) <<"Flatte_a_0(980)0: "<<physicsTree->head()->getChildValue("Flatte_a_0(980)0")
			<<"/"<<a0Res->evaluateAmp(point)*a0Res->GetNormalization();
	BOOST_LOG_TRIVIAL(info) <<"N_a_0(980)0: "<<physicsTree->head()->getChildValue("N_a_0(980)0").real()
			<<"/"<<a0Res->GetNormalization();
	BOOST_LOG_TRIVIAL(info) <<"FlatteRes_a_0(980)0: "<<physicsTree->head()->getChildValue("FlatteRes_a_0(980)0")
			<<"/"<<a0Res->evaluateAmp(point);
	BOOST_LOG_TRIVIAL(info) <<"AngD_a_0(980)0: "<<physicsTree->head()->getChildValue("AngD_a_0(980)0").real()
			<<"/"<<a0Res->evaluateWignerD(point);
	BOOST_LOG_TRIVIAL(info) <<"===========================================";

	std::shared_ptr<FitResult> resultTree = run.startFit(fitParTree);
	double finalLHTree = resultTree->getResult();
	resultTree->setTrueParameters(truePar);//set true parameters

	//======================= AMPLITUDE FIT =============================
	esti->resetInstance();
	esti = std::shared_ptr<ControlParameter>(MinLogLH::createInstance(fitAmp, inputData, toyPhspData));
	double initialLH = esti->controlParameter(fitPar);
	std::shared_ptr<Optimizer> opti(new MinuitIF(esti, fitPar));
	run.setOptimizer(opti);

	std::shared_ptr<FitResult> result= run.startFit(fitPar);
	double finalLH = result->getResult();
	result->setTrueParameters(truePar);//set true parameters

	//======================= OUTPUT =============================
	BOOST_LOG_TRIVIAL(info) <<"TREE fit result:";
	resultTree->print("P");
	BOOST_LOG_TRIVIAL(info) <<"AMPLITUDE fit result:";
	result->print("P");
	BOOST_LOG_TRIVIAL(info) <<"Comparison TREE/AMPLITUDE:";
	BOOST_LOG_TRIVIAL(info) <<"Timings[s]: "<<resultTree->getTime()<<"/"<<result->getTime();
	BOOST_LOG_TRIVIAL(info) <<"Initial likelihood: "<<initialLHTree<< "/"<<initialLH << " Deviation = "<<initialLHTree-initialLH;
	BOOST_LOG_TRIVIAL(info) <<"Final likelihood: "<<finalLHTree<< "/"<<finalLH<< " Deviation = "<<finalLHTree-finalLH;


	BOOST_LOG_TRIVIAL(info) << "FINISHED!";
	return 0;
}
