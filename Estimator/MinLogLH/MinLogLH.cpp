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

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Kinematics.hpp"

MinLogLH::MinLogLH(std::shared_ptr<Amplitude> amp_, std::shared_ptr<Amplitude> bkg_, std::shared_ptr<Data> data_,
		std::shared_ptr<Data> phspSample_,std::shared_ptr<Data> accSample_,
		unsigned int startEvent, unsigned int nEvents, double sigFrac) :
		amp(amp_), ampBkg(bkg_),data(data_), phspSample(phspSample_),accSample(accSample_), nEvts_(0), nPhsp_(0),
		nStartEvt_(startEvent), nUseEvt_(nEvents), useFunctionTree(0), signalFraction(sigFrac), accSampleEff(0), penaltyLambda(0.0)
{
	BOOST_LOG_TRIVIAL(debug)<<"MinLogLH: Constructing instance!";
	nEvts_ = data->getNEvents();
	nPhsp_ = phspSample->getNEvents();
	if(!nUseEvt_) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nEvts_)) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nPhsp_)) nUseEvt_ = nPhsp_-startEvent;
	mData = data->getMasses();
	if(data->hasWeights() && signalFraction!=1.)
		throw std::runtime_error("MinLogLH::MinLogLhbkg() data sample has weights and signal fraction is !=1. That makes no sense!");
	mPhspSample = phspSample->getMasses();
	if(accSample){
		mAccSample = accSample->getMasses();
		//we assume that the total efficiency of the sample is stored as efficiency of each event
		accSampleEff = mAccSample.eff.at(0);
		BOOST_LOG_TRIVIAL(info)<<"MinLogLH::MinLogLH() total efficiency of unbinned correction sample is set to "<<accSampleEff;
	}
	BOOST_LOG_TRIVIAL(debug)<<"MinLogLH: fraction of signal is set to "<<signalFraction<<".";
	calcSumOfWeights();
	if(signalFraction!=1 && !ampBkg)
		throw std::runtime_error("MinLogLH::MinLogLH() a signal fraction smaller 1 was set"
				" but no background description given. If you want to assume a flat background, "
				"pass a Physics//Background//FlatBackground object!");

	//reset all trees before generating new trees; saves a lot of virtual memory
	signalPhspTree = std::shared_ptr<FunctionTree>();
	signalTree_amp = std::shared_ptr<FunctionTree>();
	signalPhspTree_amp = std::shared_ptr<FunctionTree>();
	bkgPhspTree = std::shared_ptr<FunctionTree>();
	bkgTree_amp = std::shared_ptr<FunctionTree>();
	bkgPhspTree_amp = std::shared_ptr<FunctionTree>();

	calls=0;//member of ControlParameter
	return;
}

std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<Amplitude> amp_,
		std::shared_ptr<Data> data_, std::shared_ptr<Data> phspSample_,
		unsigned int startEvent, unsigned int nEvents){
	if(!instance_){
		std::shared_ptr<Data> accSample_ = std::shared_ptr<Data>();
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(
				amp_, std::shared_ptr<Amplitude>(), data_, phspSample_, accSample_, startEvent, nEvents, 1.0 ) );
		BOOST_LOG_TRIVIAL(debug)<<"MinLogLH: creating instance from amplitude and dataset!";
	}
	return instance_;
}

std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<Amplitude> amp_,
		std::shared_ptr<Data> data_, std::shared_ptr<Data> phspSample_,std::shared_ptr<Data> accSample_,
		unsigned int startEvent, unsigned int nEvents){

	if(!instance_){
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(
				amp_, std::shared_ptr<Amplitude>(), data_, phspSample_, accSample_, startEvent, nEvents, 1.) );
	}
	return instance_;
}
std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<Amplitude> amp_,std::shared_ptr<Amplitude> bkg_,
		std::shared_ptr<Data> data_, std::shared_ptr<Data> phspSample_,std::shared_ptr<Data> accSample_,
		unsigned int startEvent, unsigned int nEvents, double sigFrac){
	if(!instance_){
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(
				amp_, bkg_, data_, phspSample_, accSample_, startEvent, nEvents, sigFrac) );
	}
	return instance_;
}

void MinLogLH::setAmplitude(std::shared_ptr<Amplitude> amp_, std::shared_ptr<Data> data_,
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
	if(accSample){
		mAccSample = accSample->getMasses();
		//we assume that the total efficiency of the sample is stored as efficiency of each event
		accSampleEff = mAccSample.eff.at(0);
		BOOST_LOG_TRIVIAL(info)<<"MinLogLH::MinLogLH() total efficiency of unbinned correction "
				"sample is set to "<<accSampleEff;
	}

	signalFraction = sigFrac;
	useFunctionTree = 0;//ensure that iniLHtree is executed
	setUseFunctionTree(useFuncTr);
	if(data->hasWeights() && signalFraction!=1.)
		throw std::runtime_error("MinLogLH::MinLogLhbkg() data sample has weights and "
				"signal fraction !=1. That makes no sense!");
	calcSumOfWeights();

	if(signalFraction!=1 && !ampBkg)
		throw std::runtime_error("MinLogLH::setAmplitude() a signal fraction != 1 was set"
				" but no background description given. You should add at least a flat background!");
	if(signalFraction==1 && !ampBkg) ampBkg = std::shared_ptr<Amplitude>();

	//reset all trees before generating new trees; saves a lot of virtual memory
	signalPhspTree = std::shared_ptr<FunctionTree>();
	signalTree_amp = std::shared_ptr<FunctionTree>();
	signalPhspTree_amp = std::shared_ptr<FunctionTree>();
	bkgPhspTree = std::shared_ptr<FunctionTree>();
	bkgTree_amp = std::shared_ptr<FunctionTree>();
	bkgPhspTree_amp = std::shared_ptr<FunctionTree>();

	return;
}
void MinLogLH::setAmplitude(std::shared_ptr<Amplitude> amp_, std::shared_ptr<Amplitude> bkg_,
		std::shared_ptr<Data> data_, std::shared_ptr<Data> phspSample_, std::shared_ptr<Data> accSample_,
		unsigned int startEvent, unsigned int nEvents, bool useFuncTr, double sigFrac){
	ampBkg=bkg_;
	return setAmplitude(amp_,data_,phspSample_,accSample_,startEvent,nEvents,useFuncTr,sigFrac);
}

void MinLogLH::calcSumOfWeights(){
	sumOfWeights=0;
	for(unsigned int evt = nStartEvt_; evt<nUseEvt_+nStartEvt_; evt++){
		Event theEvent(data->getEvent(evt));
		sumOfWeights += theEvent.getWeight();
	}
	BOOST_LOG_TRIVIAL(info)<<"MinLogLH: for current data set: numEvents = "
			<<nUseEvt_<<" sumOfWeights="<<sumOfWeights<< " for current data set.";
	return;
}

void MinLogLH::setUseFunctionTree(bool t) {
	if(t==0) useFunctionTree=0;
	else iniLHtree();
}

void MinLogLH::iniLHtree(){
	//reset all trees before generating new trees; saves a lot of virtual memory
	signalPhspTree = std::shared_ptr<FunctionTree>();
	signalTree_amp = std::shared_ptr<FunctionTree>();
	signalPhspTree_amp = std::shared_ptr<FunctionTree>();
	bkgPhspTree = std::shared_ptr<FunctionTree>();
	bkgTree_amp = std::shared_ptr<FunctionTree>();
	bkgPhspTree_amp = std::shared_ptr<FunctionTree>();

	BOOST_LOG_TRIVIAL(debug) << "MinLogLH::iniLHtree() constructing the LH tree";

	if(useFunctionTree) return;
	if(!amp->hasTree()){
		throw std::runtime_error("MinLogLH::iniLHtree() amplitude has no tree");
	}
	if(ampBkg && !ampBkg->hasTree()){
		throw std::runtime_error("MinLogLH::iniLHtree() amplitude has no tree");
	}

	//----Strategies needed
	std::shared_ptr<Strategy> mmultStrat(new MultAll(ParType::MCOMPLEX));
	std::shared_ptr<Strategy> mmultDStrat(new MultAll(ParType::MDOUBLE));
	std::shared_ptr<Strategy> multiDoubleAddStrat(new AddAll(ParType::MDOUBLE));
	std::shared_ptr<Strategy> multiComplexAddStrat(new AddAll(ParType::MCOMPLEX));
	std::shared_ptr<Strategy> msqStrat(new AbsSquare(ParType::MDOUBLE));
	std::shared_ptr<Strategy> mlogStrat(new LogOf(ParType::MDOUBLE));
	std::shared_ptr<Strategy> multStrat(new MultAll(ParType::COMPLEX));
	std::shared_ptr<Strategy> multDStrat(new MultAll(ParType::DOUBLE));
	std::shared_ptr<Strategy> addStrat(new AddAll(ParType::DOUBLE));
	std::shared_ptr<Strategy> addComplexStrat(new AddAll(ParType::COMPLEX));
	std::shared_ptr<Strategy> sqStrat(new AbsSquare(ParType::DOUBLE));
	std::shared_ptr<Strategy> logStrat(new LogOf(ParType::DOUBLE));
	std::shared_ptr<Strategy> complStrat(new Complexify(ParType::COMPLEX));
	std::shared_ptr<Strategy> invStrat(new Inverse(ParType::DOUBLE));

	BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() construction normalization tree";
	//=== Signal normalization
	signalPhspTree = std::shared_ptr<FunctionTree>(new FunctionTree());
	signalPhspTree->createHead("invNormLH", invStrat);// 1/normLH
	signalPhspTree->createNode("normFactor", multDStrat, "invNormLH"); // normLH = phspVolume/N_{mc} |T_{evPHSP}|^2
	signalPhspTree->createNode("sumAmp", addStrat,"normFactor"); // sumAmp = \sum_{evPHSP} |T_{evPHSP}|^2
	std::shared_ptr<MultiDouble> eff, weightPhsp;
	signalPhspTree->createLeaf("phspVolume", Kinematics::instance()->getPhspVolume(), "normFactor");

	//Which kind of efficiency correction should be used?
	if(!accSample) {//binned
		signalPhspTree_amp = amp->GetTree(mPhspSample,mPhspSample,"_Phsp");
		signalPhspTree->createLeaf("InvNmc", 1/ ( (double) mPhspSample.sumWeight), "normFactor");
		signalPhspTree->createNode("IntensPhspEff", mmultDStrat, "sumAmp", mPhspSample.nEvents, false); //|T_{ev}|^2
		eff = std::shared_ptr<MultiDouble>( new MultiDouble("eff",mPhspSample.eff) );
		signalPhspTree->createLeaf("eff", eff, "IntensPhspEff"); //efficiency
		weightPhsp = std::shared_ptr<MultiDouble>( new MultiDouble("weightPhsp",mPhspSample.weight) );
		signalPhspTree->createLeaf("weightPhsp", weightPhsp, "IntensPhspEff"); //efficiency
		signalPhspTree->createNode("IntensPhsp", msqStrat, "IntensPhspEff", mPhspSample.nEvents, false); //|T_{ev}|^2
		BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() setting up normalization tree, "
				"using toy sample and assume that efficiency values are saved for every event!";
		//Efficiency values of toy phsp sample
		signalPhspTree->insertTree(signalPhspTree_amp, "IntensPhsp"); //Sum of resonances, at each point
	}
	else {//unbinned
		signalPhspTree->createNode("weightIntensPhsp", mmultDStrat, "sumAmp", mAccSample.nEvents, false);
		weightPhsp = std::shared_ptr<MultiDouble>( new MultiDouble("weightPhsp",mAccSample.weight) );
		signalPhspTree->createLeaf("weightPhsp", weightPhsp, "weightIntensPhsp"); //efficiency
		signalPhspTree->createNode("IntensPhsp", msqStrat, "weightIntensPhsp", mAccSample.nEvents, false); //|T_{ev}|^2
		signalPhspTree->createLeaf("InvNmc", 1/ ( (double) mAccSample.sumWeight/accSampleEff ), "normFactor");
		signalPhspTree_amp = amp->GetTree(mAccSample,mPhspSample,"_Phsp");
		BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() setting up normalization tree, "
				"using sample of accepted phsp events for efficiency correction!";
		signalPhspTree->insertTree(signalPhspTree_amp, "IntensPhsp"); //Sum of resonances, at each point
	}
	//=== Background normalization
	if(ampBkg){
		bkgPhspTree = std::shared_ptr<FunctionTree>(new FunctionTree());
		bkgPhspTree->createHead("invBkgNormLH", invStrat);// 1/normLH
		bkgPhspTree->createNode("normFactor", multDStrat, "invBkgNormLH"); // normLH = phspVolume/N_{mc} |T_{evPHSP}|^2
		bkgPhspTree->createNode("sumAmp", addStrat,"normFactor"); // sumAmp = \sum_{evPHSP} |T_{evPHSP}|^2
		bkgPhspTree->createLeaf("phspVolume", Kinematics::instance()->getPhspVolume(), "normFactor");
		if(!accSample) {//binned
			bkgPhspTree_amp = ampBkg->GetTree(mPhspSample,mPhspSample,"_Phsp");
			bkgPhspTree->createLeaf("InvNmc", 1/ ( (double) mPhspSample.sumWeight), "normFactor");
			bkgPhspTree->createNode("IntensPhspEff", mmultDStrat, "sumAmp", mPhspSample.nEvents, false); //|T_{ev}|^2
			eff = std::shared_ptr<MultiDouble>( new MultiDouble("eff",mPhspSample.eff) );
			bkgPhspTree->createLeaf("eff", eff, "IntensPhspEff"); //efficiency
			weightPhsp = std::shared_ptr<MultiDouble>( new MultiDouble("weightPhsp",mPhspSample.weight) );
			bkgPhspTree->createLeaf("weightPhsp", weightPhsp, "IntensPhspEff"); //efficiency
			bkgPhspTree->createNode("IntensPhsp", msqStrat, "IntensPhspEff", mPhspSample.nEvents, false); //|T_{ev}|^2
			BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() setting up tree for background normalization, "
					"using toy sample and assume that efficiency values are saved for every event!";
			//Efficiency values of toy phsp sample
			bkgPhspTree->insertTree(bkgPhspTree_amp, "IntensPhsp"); //Sum of resonances, at each point
		}
		else {//unbinned
			bkgPhspTree->createNode("weightIntensPhsp", mmultDStrat, "sumAmp", mAccSample.nEvents, false);
			weightPhsp = std::shared_ptr<MultiDouble>( new MultiDouble("weightPhsp",mAccSample.weight) );
			bkgPhspTree->createLeaf("weightPhsp", weightPhsp, "weightIntensPhsp"); //efficiency
			bkgPhspTree->createNode("IntensPhsp", msqStrat, "weightIntensPhsp", mAccSample.nEvents, false); //|T_{ev}|^2
			bkgPhspTree->createLeaf("InvNmc", 1/ ( (double) mAccSample.sumWeight/accSampleEff ), "normFactor");
			bkgPhspTree_amp = ampBkg->GetTree(mAccSample,mPhspSample,"_Phsp");
			BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() setting up tree for background normalization, "
					"using sample of accepted phsp events for efficiency correction!";
			bkgPhspTree->insertTree(bkgPhspTree_amp, "IntensPhsp"); //Sum of resonances, at each point
		}
	}

	BOOST_LOG_TRIVIAL(debug)<<"MinLogLH::iniLHTree() construction LH tree";
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
	signalTree_amp = amp->GetTree(mData,mPhspSample,"data");
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
	physicsTree->createNode("addBkgSignal", multiDoubleAddStrat, "Log", mData.nEvents, false);//I_{ev} = x_{ev} + (1-f_{bkg})

	//signal term
	physicsTree->createNode("normIntens", mmultDStrat, "addBkgSignal", mData.nEvents, false);// x=f_{bkg}|T|^2/norm_{LH}
	physicsTree->createLeaf("signalFrac", signalFraction, "normIntens");
	physicsTree->insertTree(signalPhspTree, "normIntens"); //provides 1/normLH
	physicsTree->createNode("Intens", msqStrat, "normIntens", mData.nEvents, false);
	physicsTree->insertTree(signalTree_amp,"Intens");
	//background term
	if(ampBkg){
		physicsTree->createNode("normBkg", mmultDStrat, "addBkgSignal", mData.nEvents, false);// x=f_{bkg}|T|^2/norm_{LH}
		physicsTree->createLeaf("OneMinusBkgFrac", (1-signalFraction), "normBkg");
		bkgTree_amp= ampBkg->GetTree(mData,mPhspSample,"data");
		physicsTree->insertTree(bkgPhspTree, "normBkg"); //provides 1/normLH
		physicsTree->createNode("IntensBkg", msqStrat, "normBkg", mData.nEvents, false);
		physicsTree->insertTree(bkgTree_amp,"IntensBkg");
	}
	physicsTree->recalculate();
	//	std::string treeString = physicsTree->head()->to_str(10);
	//	BOOST_LOG_TRIVIAL(debug) << std::endl << treeString;
	if(!physicsTree->sanityCheck()) {
		throw std::runtime_error("MinLogLH::iniLHtree() tree has structural problems. Sanity check not passed!");
	}
	BOOST_LOG_TRIVIAL(debug) <<"MinLogLH::iniLHtree() construction of LH tree finished!";
	useFunctionTree=1;
	return;
}

double MinLogLH::calcPenalty(){
	if(penaltyLambda<=0) return 0; //penalty term disabled
	double magSum = 0;
	auto it = amp->GetResonanceItrFirst();
	for(; it != amp->GetResonanceItrLast(); ++it){
		magSum += std::fabs( (*it)->GetMagnitude()) *
				std::sqrt( (*it)->GetIntegral() );
	}
	return (penaltyLambda*magSum);
}

double MinLogLH::controlParameter(ParameterList& minPar){
	amp->setParameterList(minPar); //setting new parameters

	double lh=0;
	if(!useFunctionTree){
		//Calculate normalization
		double vol = Kinematics::instance()->getPhspVolume();
		double norm=0, normSq=0, normBkg=0;
		std::shared_ptr<Data> sam;
		if(accSample) sam=accSample;
		else sam=phspSample;
		int sam_size = sam->getNEvents();
		for(unsigned int phsp=0; phsp<sam_size; phsp++){ //loop over phspSample
			Event theEvent(sam->getEvent(phsp));
			dataPoint point(theEvent);
			double intens = 0, intensBkg = 0;
			ParameterList intensL = amp->intensity(point);
			intens = intensL.GetDoubleParameter(0)->GetValue();
			if(intens>0) {
				norm+=intens;
				normSq+=intens*intens;
			}

			if(ampBkg){
				ParameterList intensB = ampBkg->intensity(point);
				intensBkg = intensB.GetDoubleParameter(0)->GetValue();
			}else{
				intensBkg = 0;
			}
			if(intensBkg>0) normBkg+=intensBkg;
		}
		normBkg = normBkg * vol/sam_size;
		if(normBkg==0) normBkg=1;
		norm = norm * vol/sam_size;
		//Approximate error of integration, see (Numerical Recipes Vol3, p398, Eq. 7.7.1)
		double normError = sqrt( ( vol*vol*normSq/sam_size - norm*norm) / sam_size);

		//Use internal amplitude integration - no unbinned efficiency correction possible
		//if(ampBkg) normBkg=ampBkg->normalization();
		//else normBkg=1;
		//norm = amp->normalization();

		if(norm==0) norm=1;
		//Calculate \Sum_{ev} log()
		double sumLog=0;
		for(unsigned int evt = nStartEvt_; evt<nUseEvt_+nStartEvt_; evt++){//loop over data sample
			Event theEvent(data->getEvent(evt));
			dataPoint point(theEvent);
			double intens = 0, intensBkg = 1;
			/*Get intensity of amplitude at data point. Efficiency is a constant term in LH and
			 *therefore we drop it here to be consistent with the tree */
			ParameterList intensS = amp->intensityNoEff(point);
			intens = intensS.GetDoubleParameter(0)->GetValue();
			if(ampBkg){
				ParameterList intensB = ampBkg->intensityNoEff(point);
				intensBkg = intensB.GetDoubleParameter(0)->GetValue();
			}else{
				intensBkg = 0;
			}
			if(intens>0) sumLog += std::log(
					signalFraction*intens/norm+(1-signalFraction)*intensBkg/normBkg )*theEvent.getWeight();
		}
		lh = (-1)*((double)nUseEvt_)/sumOfWeights*sumLog;
	} else {
		physicsTree->recalculate();
		std::shared_ptr<DoubleParameter> logLH = std::dynamic_pointer_cast<DoubleParameter>(
				physicsTree->head()->getValue() );
		lh = logLH->GetValue();
	}
	lh += calcPenalty();
	calls++;
	return lh; //return -logLH
}

void MinLogLH::setPenaltyScale(double sc) {
	if(sc < 0){
		BOOST_LOG_TRIVIAL(info) << "MinLogLH::setPenaltyScale | penalty scale cannot be negative!";
		return;
	}
	BOOST_LOG_TRIVIAL(info) << "MinLogLH::setPenaltyScale | Setting scale of penalty term to "<<sc;
	penaltyLambda = sc;
}
