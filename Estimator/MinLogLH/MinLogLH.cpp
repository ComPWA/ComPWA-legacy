//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Kinematics.hpp"

void MinLogLH::calcSumOfWeights(){
	sumOfWeights=0;
	if(pDIF_) {//if we have a data sample sum up all weights
		for(unsigned int evt = nStartEvt_; evt<nUseEvt_+nStartEvt_; evt++){
			Event theEvent(pDIF_->getEvent(evt));
			sumOfWeights += theEvent.getWeight();
		}
	} else if(pEvtTree_){
		// need to calc weights from tree node
		std::vector<double> weights = pEvtTree_->head()->getChildMultiDoubleValue("weight");
		for(unsigned int i=0; i<weights.size(); i++) sumOfWeights+=weights.at(i);
//		sumOfWeights = nUseEvt_;
	} else {
		BOOST_LOG_TRIVIAL(error)<<"MinLogLH: not data sample available. Can't calculate sumOfWeights. Using number of events!";
		sumOfWeights = nUseEvt_;
	}
	BOOST_LOG_TRIVIAL(info)<<"MinLogLH: for current data set: numEvents = "<<nUseEvt_<<" sumOfWeights="<<sumOfWeights<< " for current data set.";
	return;
}

MinLogLH::MinLogLH(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF, unsigned int startEvent, unsigned int nEvents)
: pPIF_(inPIF), pDIF_(inDIF), nEvts_(0), nPhsp_(0), nStartEvt_(startEvent), nUseEvt_(nEvents){
	//	phspVolume = Kinematics::instance()->getPhspVolume();
	nEvts_ = pDIF_->getNEvents();
	if( !(startEvent+nUseEvt_<=nEvts_) ) nUseEvt_ = nEvts_-startEvent;
	calcSumOfWeights();
}

MinLogLH::MinLogLH(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF,
		std::shared_ptr<Data> inPHSP, unsigned int startEvent, unsigned int nEvents) :
						pPIF_(inPIF), pDIF_(inDIF), pPHSP_(inPHSP), nEvts_(0), nPhsp_(0),
						nStartEvt_(startEvent), nUseEvt_(nEvents){
	//	phspVolume = Kinematics::instance()->getPhspVolume();
	nEvts_ = pDIF_->getNEvents();
	nPhsp_ = inPHSP->getNEvents();
	if(!nUseEvt_) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nEvts_)) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nPhsp_)) nUseEvt_ = nPhsp_-startEvent;
	calcSumOfWeights();
}

MinLogLH::MinLogLH(std::shared_ptr<FunctionTree> inEvtTree, unsigned int inNEvts)
: pEvtTree_(inEvtTree), nEvts_(inNEvts), nPhsp_(0), nStartEvt_(0), nUseEvt_(inNEvts){
	//	phspVolume = Kinematics::instance()->getPhspVolume();
	calcSumOfWeights();
}

MinLogLH::MinLogLH(std::shared_ptr<FunctionTree> inEvtTree, std::shared_ptr<FunctionTree> inPhspTree, unsigned int inNEvts)
: pEvtTree_(inEvtTree), pPhspTree_(inPhspTree), nEvts_(inNEvts), nPhsp_(0), nStartEvt_(0), nUseEvt_(inNEvts){
	calcSumOfWeights();

}

void MinLogLH::setTree(std::shared_ptr<FunctionTree> inEvtTree, unsigned int inNEvts){
	pEvtTree_=inEvtTree;
	nEvts_=inNEvts;
	pPhspTree_=std::shared_ptr<FunctionTree>();
	nPhsp_=0;
	calcSumOfWeights();
	return;
}
void MinLogLH::setTree(std::shared_ptr<FunctionTree> inEvtTree, std::shared_ptr<FunctionTree> inPhspTree, unsigned int inNEvts){
	pEvtTree_=inEvtTree;
	nEvts_=inNEvts;
	nUseEvt_=inNEvts;
	pPhspTree_=inPhspTree;
	calcSumOfWeights();
	nPhsp_=0;
}
void MinLogLH::setAmplitude(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF,
		std::shared_ptr<Data> inPHSP, unsigned int startEvent, unsigned int nEvents){
	pPIF_ = inPIF;
	pDIF_=inDIF;
	pPHSP_=inPHSP;
	//nEvts_=0;
	//nPhsp_=0;
	nStartEvt_= startEvent;
	nUseEvt_= nEvents;
	//	phspVolume = Kinematics::instance()->getPhspVolume();
	nEvts_ = pDIF_->getNEvents();
	nPhsp_ = inPHSP->getNEvents();
	if(!nUseEvt_) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nEvts_)) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nPhsp_)) nUseEvt_ = nPhsp_-startEvent;
	calcSumOfWeights();
	return;
}

std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<Amplitude> inPIF,
		std::shared_ptr<Data> inDIF, unsigned int startEvent, unsigned int nEvents){
	if(!instance_){
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(inPIF, inDIF, startEvent, nEvents));
		BOOST_LOG_TRIVIAL(debug)<<"MinLogLH: creating instance from amplitude and dataset!";
	}

	return instance_;
}

std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<Amplitude> inPIF,
		std::shared_ptr<Data> inDIF, std::shared_ptr<Data> inPHSP,
		unsigned int startEvent, unsigned int nEvents){
	if(!instance_){
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(inPIF, inDIF, inPHSP, startEvent, nEvents));
		BOOST_LOG_TRIVIAL(debug)<<"MinLogLH: creating instance from amplitude and dataset!";
	}

	return instance_;
}

std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<FunctionTree> inEvtTree,
		unsigned int inNEvts){
	if(!instance_){
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(inEvtTree, inNEvts));
		BOOST_LOG_TRIVIAL(debug)<<"MinLogLH: creating instance from single tree!";
	}

	return instance_;
}

std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<FunctionTree> inEvtTree,
		std::shared_ptr<FunctionTree> inPhspTree, unsigned int inNEvts){
	if(!instance_){
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(inEvtTree, inPhspTree, inNEvts));
		BOOST_LOG_TRIVIAL(debug)<<"MinLogLH: creating instance from two trees!";
	}
	return instance_;
}

MinLogLH::~MinLogLH(){
	//delete instance_;
}

double MinLogLH::controlParameter(ParameterList& minPar){

	if(pPIF_) pPIF_->setParameterList(minPar); //setting new parameters

	//Calculate normalization
	double norm = 0;
	if(pPHSP_&& pPIF_){//norm by phasespace monte-carlo
		//		for(unsigned int phsp=nStartEvt_; phsp<nUseEvt_+nStartEvt_; phsp++){//TODO: needs review
		for(unsigned int phsp=0; phsp<nPhsp_; phsp++){
			Event theEvent(pPHSP_->getEvent(phsp));
			if(theEvent.getNParticles()!=3) continue;
			dataPoint point(theEvent);
			double intens = 0;
			ParameterList intensL = pPIF_->intensity(point);
			intens = intensL.GetDoubleParameter(0)->GetValue();
			if(intens>0) norm+=intens;
		}
	}else if(pPhspTree_){//normalization by tree
		pPhspTree_->recalculate();
		std::shared_ptr<DoubleParameter> intensL = std::dynamic_pointer_cast<DoubleParameter>(pPhspTree_->head()->getValue());
		norm = intensL->GetValue();
	}else if(!pPHSP_ && pPIF_){//normalization by numerical integration
		norm=pPIF_->integral();
	}else{
		BOOST_LOG_TRIVIAL(error)<< "MinLogLH::controlParameter() : no tree and no amplitude given. Can't calculate Normalization!";
		//TODO: Exception
		return 0;
	}

	//Calculate likelihood
	double lh=0; //calculate LH:
	if(pDIF_ && pPIF_){//amplitude and datasample
		double f = 1.0;
		for(unsigned int evt = nStartEvt_; evt<nUseEvt_+nStartEvt_; evt++){
			Event theEvent(pDIF_->getEvent(evt));
			dataPoint point(theEvent);
			double intens = 0;
			/*Get intensity of amplitude at data point. Efficiency is a constant term in LH and
			 *therefore we drop it here to be consistent with the tree */
			ParameterList intensL = pPIF_->intensityNoEff(point);
			intens = intensL.GetDoubleParameter(0)->GetValue();
			if(intens>0) lh += std::log(intens)*theEvent.getWeight();
//			if(intens>0) lh += std::log( f*intens/norm+(1-f) )*theEvent.getWeight();
		}
//		lh = sumOfWeights*std::log(norm) - lh ;//other factors are constant and drop in deviation, so we can ignore them
//		lh = std::log(norm) - lh/sumOfWeights ;//other factors are constant and drop in deviation, so we can ignore them
		lh = (-1)*nUseEvt_/sumOfWeights*lh ;//other factors are constant and drop in deviation, so we can ignore them
	}else if(pEvtTree_){//tree
		pEvtTree_->recalculate();
		std::shared_ptr<DoubleParameter> intensL = std::dynamic_pointer_cast<DoubleParameter>(pEvtTree_->head()->getValue());
		lh = intensL->GetValue();
//		lh = sumOfWeights*std::log(norm) - lh ;//other factors are constant and drop in deviation, so we can ignore them
//		lh = std::log(norm)*nUseEvt_ - lh*nUseEvt_/sumOfWeights ;//other factors are constant and drop in deviation, so we can ignore them
	}else{
		BOOST_LOG_TRIVIAL(error)<< "MinLogLH::controlParameter() : no tree and no amplitude given. Can't calculate LH!";
		//TODO: Exception
		return 0;
	}

//	BOOST_LOG_TRIVIAL(debug) << "Data Term: " << lh << "\t Phsp Term (wo log): " << norm;
//	lh = nUseEvt_*std::log(norm) - lh ;//other factors are constant and drop in deviation, so we can ignore them
	lh = sumOfWeights*std::log(norm) - lh ;//other factors are constant and drop in deviation, so we can ignore them

	return lh; //return -logLH
}
