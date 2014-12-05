//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//	   Peter Weidenkaff - Weights and background fractions
//-------------------------------------------------------------------------------
#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <exception>

#include "Estimator/MinLogLHbkg/MinLogLHbkg.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Kinematics.hpp"

MinLogLHbkg::MinLogLHbkg(std::shared_ptr<Amplitude> amp_, std::shared_ptr<Data> data_,
		std::shared_ptr<Data> phspSample_,std::shared_ptr<Data> accSample_,
		unsigned int startEvent, unsigned int nEvents, double sigFrac) :
		amp(amp_), data(data_), phspSample(phspSample_),accSample(accSample_), nEvts_(0), nPhsp_(0),
		nStartEvt_(startEvent), nUseEvt_(nEvents), useFunctionTree(0), signalFraction(sigFrac){
	BOOST_LOG_TRIVIAL(debug)<<"MinLogLHbkg: Constructing instance!";
	nEvts_ = data->getNEvents();
	nPhsp_ = phspSample->getNEvents();
	if(!nUseEvt_) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nEvts_)) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nPhsp_)) nUseEvt_ = nPhsp_-startEvent;
	mData = data->getMasses();
	if(data->hasWeights() && signalFraction!=1.)
		throw std::runtime_error("MinLogLHbkg::MinLogLhbkg() data sample has weights and signal fraction is !=1. That makes no sense!");
	mPhspSample = phspSample->getMasses();
	if(accSample) mAccSample = accSample->getMasses();
	BOOST_LOG_TRIVIAL(debug)<<"MinLogLHbkg: fraction of signal is set to "<<signalFraction<<".";
	calcSumOfWeights();
}

std::shared_ptr<ControlParameter> MinLogLHbkg::createInstance(std::shared_ptr<Amplitude> amp_,
		std::shared_ptr<Data> data_, std::shared_ptr<Data> phspSample_,
		unsigned int startEvent, unsigned int nEvents, double sigFrac){
	if(!instance_){
		std::shared_ptr<Data> accSample_ = std::shared_ptr<Data>();
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLHbkg(
				amp_, data_, phspSample_, accSample_, startEvent, nEvents, sigFrac) );
		BOOST_LOG_TRIVIAL(debug)<<"MinLogLHbkg: creating instance from amplitude and dataset!";
	}
	return instance_;
}

std::shared_ptr<ControlParameter> MinLogLHbkg::createInstance(std::shared_ptr<Amplitude> amp_,
		std::shared_ptr<Data> data_, std::shared_ptr<Data> phspSample_,std::shared_ptr<Data> accSample_,
		unsigned int startEvent, unsigned int nEvents, double sigFrac){
	if(!instance_){
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLHbkg(
				amp_, data_, phspSample_, accSample_, startEvent, nEvents, sigFrac) );
	}
	return instance_;
}

void MinLogLHbkg::setAmplitude(std::shared_ptr<Amplitude> amp_, std::shared_ptr<Data> data_,
		std::shared_ptr<Data> phspSample_, std::shared_ptr<Data> accSample_,
		unsigned int startEvent, unsigned int nEvents, bool useFuncTr, double sigFrac){
	amp = amp_;
	data = data_;
	phspSample = phspSample_;
	accSample = accSample_;

	nStartEvt_= startEvent;
	nUseEvt_= nEvents;
	nEvts_ = data->getNEvents();
	nPhsp_ = phspSample->getNEvents();
	if(!nUseEvt_) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nEvts_)) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nPhsp_)) nUseEvt_ = nPhsp_-startEvent;
	mData = data->getMasses();
	mPhspSample = phspSample->getMasses();
	if(accSample) mAccSample = accSample->getMasses();
	signalFraction = sigFrac;
	useFunctionTree = 0;//ensure that iniLHtree is executed
	setUseFunctionTree(useFuncTr);
	if(data->hasWeights() && signalFraction!=1.)
		throw std::runtime_error("MinLogLHbkg::MinLogLhbkg() data sample has weights and signal fraction !=1. That makes no sense!");
	calcSumOfWeights();
	return;
}

void MinLogLHbkg::calcSumOfWeights(){
	sumOfWeights=0;
	for(unsigned int evt = nStartEvt_; evt<nUseEvt_+nStartEvt_; evt++){
		Event theEvent(data->getEvent(evt));
		sumOfWeights += theEvent.getWeight();
	}
	BOOST_LOG_TRIVIAL(info)<<"MinLogLHbkg: for current data set: numEvents = "<<nUseEvt_<<" sumOfWeights="<<sumOfWeights<< " for current data set.";
	return;
}

void MinLogLHbkg::setUseFunctionTree(bool t) {
	if(t==0) useFunctionTree=0;
	else iniLHtree();
}

void MinLogLHbkg::iniLHtree(){
	BOOST_LOG_TRIVIAL(debug) << "MinLogLHbkg::iniLHtree() constructing the LH tree";

	if(useFunctionTree) return;
	if(!amp->hasTree()){
		throw std::runtime_error("MinLogLHbkg::iniLHtree() amplitude has no tree");
	}

	//----Strategies needed
	std::shared_ptr<MultAll> mmultStrat = std::shared_ptr<MultAll>(new MultAll(ParType::MCOMPLEX));
	std::shared_ptr<MultAll> mmultDStrat = std::shared_ptr<MultAll>(new MultAll(ParType::MDOUBLE));
	std::shared_ptr<AddAll> multiDoubleAddStrat = std::shared_ptr<AddAll>(new AddAll(ParType::MDOUBLE));
	std::shared_ptr<AddAll> multiComplexAddStrat = std::shared_ptr<AddAll>(new AddAll(ParType::MCOMPLEX));
	std::shared_ptr<AbsSquare> msqStrat = std::shared_ptr<AbsSquare>(new AbsSquare(ParType::MDOUBLE));
	std::shared_ptr<LogOf> mlogStrat = std::shared_ptr<LogOf>(new LogOf(ParType::MDOUBLE));
	std::shared_ptr<MultAll> multStrat = std::shared_ptr<MultAll>(new MultAll(ParType::COMPLEX));
	std::shared_ptr<MultAll> multDStrat = std::shared_ptr<MultAll>(new MultAll(ParType::DOUBLE));
	std::shared_ptr<AddAll> addStrat = std::shared_ptr<AddAll>(new AddAll(ParType::DOUBLE));
	std::shared_ptr<AddAll> addComplexStrat = std::shared_ptr<AddAll>(new AddAll(ParType::COMPLEX));
	std::shared_ptr<AbsSquare> sqStrat = std::shared_ptr<AbsSquare>(new AbsSquare(ParType::DOUBLE));
	std::shared_ptr<LogOf> logStrat = std::shared_ptr<LogOf>(new LogOf(ParType::DOUBLE));
	std::shared_ptr<Complexify> complStrat = std::shared_ptr<Complexify>(new Complexify(ParType::COMPLEX));
	std::shared_ptr<Inverse> invStrat = std::shared_ptr<Inverse>(new Inverse(ParType::DOUBLE));
	std::shared_ptr<SquareRoot> sqRootStrat = std::shared_ptr<SquareRoot>(new SquareRoot(ParType::DOUBLE));

	BOOST_LOG_TRIVIAL(debug)<<"MinLogLHbkg::iniLHTree() construction normalization tree";
	//====== PHSP tree for LH normalization
	signalPhspTree = std::shared_ptr<FunctionTree>(new FunctionTree());
	signalPhspTree->createHead("invNormLH", invStrat);// 1/normLH
	signalPhspTree->createNode("normFactor", multDStrat, "invNormLH"); // normLH = phspVolume/N_{mc} |T_{evPHSP}|^2
	signalPhspTree->createLeaf("phspVolume", Kinematics::instance()->getPhspVolume(), "normFactor");
	signalPhspTree->createNode("InvNmc", invStrat, "normFactor");
	signalPhspTree->createLeaf("Nmc", mPhspSample.nEvents, "InvNmc");
	signalPhspTree->createNode("sumAmp", addStrat,"normFactor"); // sumAmp = \sum_{evPHSP} |T_{evPHSP}|^2
	std::shared_ptr<MultiDouble> eff;

	//Which kind of efficiency correction should be used?
	if(!accSample) {//binned
		signalPhspTree_amp = amp->getAmpTree(mPhspSample,mPhspSample,"_Phsp");
		signalPhspTree->createNode("IntensPhsp", msqStrat, "sumAmp", mPhspSample.nEvents, false); //|T_{ev}|^2
		BOOST_LOG_TRIVIAL(debug)<<"MinLogLHbkg::iniLHTree() setting up normalization tree, "
				"using toy sample and assume that efficiency values are saved for every event!";
		//Efficiency values of toy phsp sample
		eff = std::shared_ptr<MultiDouble>( new MultiDouble("eff",mPhspSample.eff) ); //only needed for opt == "norm"
		signalPhspTree->createNode("AmplitudeEff", mmultStrat, "IntensPhsp", mPhspSample.nEvents, false); //Sum of resonances, at each point
		signalPhspTree->createLeaf("eff", eff, "AmplitudeEff"); //efficiency
		signalPhspTree->insertTree(signalPhspTree_amp, "AmplitudeEff"); //Sum of resonances, at each point
	}
	else {//unbinned
		signalPhspTree->createNode("IntensPhsp", msqStrat, "sumAmp", mAccSample.nEvents, false); //|T_{ev}|^2
		signalPhspTree_amp = amp->getAmpTree(mAccSample,mPhspSample,"_Phsp");
		BOOST_LOG_TRIVIAL(debug)<<"MinLogLHbkg::iniLHTree() setting up normalization tree, "
				"using sample of accepted phsp events for efficiency correction!";
		signalPhspTree->insertTree(signalPhspTree_amp, "IntensPhsp"); //Sum of resonances, at each point
	}

	BOOST_LOG_TRIVIAL(debug)<<"MinLogLHbkg::iniLHTree() construction LH tree";
	/* CONSTRUCTION OF THE LIKELIHOOD:
	 * We denote the coherent sum over all resonances with T:
	 * 		T := \sum_{i,j} c_i c_j^*A_iA_j^*
	 * The negative log LH is given by:
	 * 		-log L = - N/(\sum_{ev} w_{ev}) \sum_{ev} w_{ev} \log{f_{bkg} \frac{|T|^2}{\int_{DP} |T|^2} + (1-f_{bkg})}
	 * The sum over all weights is necessary to normalize the weights to one. Otherwise the error
	 * estimate is incorrect. The LH normalization is norm_{LH} = \int_{DP} |T|^2.
	 * This formulation includes event weights as well as a flat background description. f_{bkg} is
	 * the fraction of background in the sample. Using both is of course non-sense. Set weights to
	 * one OR f_{bkg} to zero.
	 */
	physicsTree = std::shared_ptr<FunctionTree>(new FunctionTree());
	/* Setup basic tree
	 * head node is 'Amplitude' which contains the complex amplitude values for each event in sample
	 */
	signalTree_amp = amp->getAmpTree(mData,mPhspSample,"");
	//------------Setup Tree Pars---------------------
	std::shared_ptr<MultiDouble> weight = std::shared_ptr<MultiDouble>( new MultiDouble("weight",mData.weight) );

	physicsTree->createHead("LH", multDStrat); //-log L = (-1)*N/(\sum_{ev} w_{ev}) \sum_{ev} ...
	physicsTree->createLeaf("minusOne", -1 ,"LH");
	physicsTree->createLeaf("nEvents", mData.nEvents ,"LH");
	physicsTree->createNode("invSumWeights", invStrat,"LH"); // 1/\sum_{ev} w_{ev}
	physicsTree->createNode("sumEvents", addStrat, "LH"); // \sum_{ev} w_{ev} * log( I_{ev} )
	physicsTree->createNode("sumWeights", addStrat, "invSumWeights"); // \sum_{ev} w_{ev}
	physicsTree->createLeaf("weight", weight, "sumWeights");
	physicsTree->createNode("weightLog", mmultDStrat, "sumEvents", mData.nEvents, false); //w_{ev} * log( I_{ev} )
	physicsTree->createLeaf("weight", weight, "weightLog");
	physicsTree->createNode("Log", mlogStrat, "weightLog", mData.nEvents, false); // log(I_{ev})
	physicsTree->createNode("addBkg", multiDoubleAddStrat, "Log", mData.nEvents, false);//I_{ev} = x_{ev} + (1-f_{bkg})
	physicsTree->createLeaf("OneMinusBkgFrac", (1-signalFraction), "addBkg");
	physicsTree->createNode("normIntens", mmultDStrat, "addBkg", mData.nEvents, false);// x=f_{bkg}|T|^2/norm_{LH}
	physicsTree->insertTree(signalPhspTree, "normIntens"); //provides 1/normLH
	physicsTree->createLeaf("signalFrac", signalFraction, "normIntens");
	physicsTree->createNode("Intens", msqStrat, "normIntens", mData.nEvents, false);
	physicsTree->insertTree(signalTree_amp,"Intens");

	BOOST_LOG_TRIVIAL(debug) << std::endl << physicsTree;
	if(!physicsTree->sanityCheck()) {
		throw std::runtime_error("MinLogLHbkg::iniLHtree() tree has structural problems. Sanity check not passed!");
	}
	BOOST_LOG_TRIVIAL(debug) <<"MinLogLHbkg::iniLHtree() construction of LH tree finished!";
	useFunctionTree=1;
	return;
}

double MinLogLHbkg::controlParameter(ParameterList& minPar){
	amp->setParameterList(minPar); //setting new parameters

	double lh=0, norm=0;
	if(!useFunctionTree){
		//Calculate normalization
		for(unsigned int phsp=0; phsp<nPhsp_; phsp++){ //loop over phspSample
			Event theEvent(phspSample->getEvent(phsp));
			if(theEvent.getNParticles()!=3) continue;
			dataPoint point(theEvent);
			double intens = 0;
			ParameterList intensL = amp->intensity(point);
			intens = intensL.GetDoubleParameter(0)->GetValue();
			if(intens>0) norm+=intens;
		}
		//Calculate -logLH
		for(unsigned int evt = nStartEvt_; evt<nUseEvt_+nStartEvt_; evt++){//loop over data sample
			Event theEvent(data->getEvent(evt));
			dataPoint point(theEvent);
			double intens = 0;
			/*Get intensity of amplitude at data point. Efficiency is a constant term in LH and
			 *therefore we drop it here to be consistent with the tree */
			ParameterList intensL = amp->intensityNoEff(point);
			intens = intensL.GetDoubleParameter(0)->GetValue();
			//			if(intens>0) lh += std::log(intens)*theEvent.getWeight();
			if(intens>0) lh += std::log( signalFraction*intens/norm+(1-signalFraction) )*theEvent.getWeight();
		}
		lh = (-1)*nUseEvt_/sumOfWeights*lh ;//other factors are constant and drop in deviation, so we can ignore them
	} else {
		physicsTree->recalculate();
		std::shared_ptr<DoubleParameter> logLH = std::dynamic_pointer_cast<DoubleParameter>(
				physicsTree->head()->getValue() );
		lh = logLH->GetValue();
	}
	//BOOST_LOG_TRIVIAL(debug) << "-logLH: " << lh ;

	return lh; //return -logLH
}
