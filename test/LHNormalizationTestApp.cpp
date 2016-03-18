//-------------------------------------------------------------------------------
// Copyright (c) 2015 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff
//-------------------------------------------------------------------------------
/** Test likelihood normalization of MinLogLH
 * @file LHNormalizationTestApp.cpp
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

// ComPWA header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/RunManager.hpp"
#include "Core/Efficiency.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/TableFormater.hpp"
#include "Core/AbsParameter.hpp"
#include "Core/Logging.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Physics/DPKinematics/RootEfficiency.cpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Physics/DPKinematics/RootGenerator.hpp"

using namespace std;

/************************************************************************************************/
void randomStartValues(ParameterList& fitPar){
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
//	std::cout<<"Randomizing parameter list. New list:"<<fitPar<<std::endl;
	return;
}
/**
 * The main function.
 */
int main(int argc, char **argv){
	int seed = 3041; //default seed

	//number of calls for numeric integration and number of events for phsp integration
	unsigned int num = 1000000;
	Logging log("log-compareTreeAmp.txt",boost::log::trivial::debug); //initialize logging
	//initialize kinematics of decay
	DalitzKinematics::createInstance("D0","K_S0","K-","K+");//setup kinematics
	//initialize random generator
	std::shared_ptr<Generator> gen = std::shared_ptr<Generator>(new RootGenerator(seed));

	std::cout<<std::setprecision(8)<<std::endl;
	RunManager run;
	run.setGenerator(gen);
	std::shared_ptr<Efficiency> eff(new UnitEfficiency());

	std::shared_ptr<Data> toyPhspData(new RootReader());//empty phsp sample
	run.setPhspSample(toyPhspData);
	if( !toyPhspData->getNEvents() ) {
		run.generatePhsp(num);
		num=toyPhspData->getNEvents();
	}

	std::shared_ptr<Amplitude> unitAmp(new UnitAmp());
	ParameterList list;
	std::shared_ptr<ControlParameter> esti;
	//Use unit amplitude
	//	esti = std::shared_ptr<ControlParameter>(MinLogLH::createInstance(unitAmp,
	//		toyPhspData, toyPhspData));


	//Use example model with 3 resonances
	std::string trueModelFile = "test/CompareTreeAmp-model.xml";
	boost::property_tree::ptree pt;
	read_xml(trueModelFile, pt, boost::property_tree::xml_parser::trim_whitespace);
	auto fitAmpPtr = new AmpSumIntensity("amp",normStyle::none, eff, num);
	fitAmpPtr->Configure(pt);
	std::shared_ptr<Amplitude> trueAmp( fitAmpPtr );
	trueAmp->FillParameterList(list);
	esti = std::shared_ptr<ControlParameter>(MinLogLH::createInstance(trueAmp,toyPhspData , toyPhspData));

	MinLogLH* minLog = dynamic_cast<MinLogLH*>(&*(esti->Instance()));
	minLog->setUseFunctionTree(1);
	std::shared_ptr<FunctionTree> physicsTree = minLog->getTree();
	BOOST_LOG_TRIVIAL(debug) << physicsTree->head()->to_str(10);
	double initialLHTree = esti->controlParameter(list);
	std::shared_ptr<Optimizer> optiTree(new MinuitIF(esti, list));
	run.setOptimizer(optiTree);

	for(int i=0; i<100; i++){
		if(i>0)	randomStartValues(list);
		trueAmp->UpdateParameters(list);
		esti->controlParameter(list);
		std::shared_ptr<MultiDouble> in = std::dynamic_pointer_cast<MultiDouble>(
				physicsTree->head()->getChildValue("Intens"));
		double sumNormAmp = physicsTree->head()->getChildSingleValue("sumAmp").real();
		double norm = physicsTree->head()->getChildSingleValue("normFactor").real();
		double norm_expect = Kinematics::instance()->GetPhspVolume()/num*sumNormAmp;
		double sumLog = 0;
		for(int i=0; i<in->GetNValues();i++)
			sumLog+=std::log(in->getNodeValue(i).real());
		double lh_expect = sumLog
				-num*std::log(Kinematics::instance()->GetPhspVolume()/num)
		-num*std::log(sumNormAmp);
		lh_expect*=-1;
		//	std::cout<<sumLog<<" "<<n*std::log(Kinematics::instance()->getPhspVolume()/n)
		//	<<" "<<n*std::log(sumNormAmp)<<std::endl;
		std::shared_ptr<DoubleParameter> intens = std::dynamic_pointer_cast<DoubleParameter>(
				physicsTree->head()->getValue());

		BOOST_LOG_TRIVIAL(info) <<"Deviation LH(value - expectation): "
				<<intens->GetValue()-lh_expect
				<<"    Deviation norm: " <<(norm - norm_expect);
	}

	BOOST_LOG_TRIVIAL(info) << "FINISHED!";
	return 0;
}
