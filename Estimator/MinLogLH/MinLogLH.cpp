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
//#include "Physics/DPKinematics/DataPoint.hpp"

MinLogLH::MinLogLH(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF)
: pPIF_(inPIF), pDIF_(inDIF), nEvts_(0), nPhsp_(0){
phspVolume = Kinematics::instance()->getPhspVolume();
nEvts_ = pDIF_->getNEvents();

}

MinLogLH::MinLogLH(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF, std::shared_ptr<Data> inPHSP)
: pPIF_(inPIF), pDIF_(inDIF), pPHSP_(inPHSP), nEvts_(0), nPhsp_(0){
phspVolume = Kinematics::instance()->getPhspVolume();
nEvts_ = pDIF_->getNEvents();
nPhsp_ = inPHSP->getNEvents();
}

MinLogLH::MinLogLH(std::shared_ptr<FunctionTree> inEvtTree, unsigned int inNEvts)
: pEvtTree_(inEvtTree), nEvts_(inNEvts), nPhsp_(0){
phspVolume = Kinematics::instance()->getPhspVolume();

}

MinLogLH::MinLogLH(std::shared_ptr<FunctionTree> inEvtTree, std::shared_ptr<FunctionTree> inPhspTree, unsigned int inNEvts, unsigned int inNPhsp)
: pEvtTree_(inEvtTree), pPhspTree_(inPhspTree), nEvts_(inNEvts), nPhsp_(inNPhsp){
phspVolume = Kinematics::instance()->getPhspVolume();

}

std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF){
	if(!instance_)
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(inPIF, inDIF));

	return instance_;
}

std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<Amplitude> inPIF, std::shared_ptr<Data> inDIF, std::shared_ptr<Data> inPHSP){
	if(!instance_)
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(inPIF, inDIF, inPHSP ));

	return instance_;
}

std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<FunctionTree> inEvtTree, unsigned int inNEvts){
	if(!instance_)
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(inEvtTree, inNEvts));

	return instance_;
}

std::shared_ptr<ControlParameter> MinLogLH::createInstance(std::shared_ptr<FunctionTree> inEvtTree, std::shared_ptr<FunctionTree> inPhspTree, unsigned int inNEvts, unsigned int inNPhsp){
	if(!instance_)
		instance_ = std::shared_ptr<ControlParameter>(new MinLogLH(inEvtTree, inPhspTree, inNEvts, inNPhsp));

	return instance_;
}

MinLogLH::~MinLogLH(){
	//delete instance_;
}

double MinLogLH::controlParameter(ParameterList& minPar){
	/*unsigned int nEvents = pDIF_->getNEvents();
	unsigned int nPHSPEvts=0;
	if(pPHSP_) nPHSPEvts = pPHSP_->getNEvents();
	unsigned int nParts = ((Event)pDIF_->getEvent(0)).getNParticles();*/

	//check if able to handle this many particles
	/*if(nParts<2 || nParts>3){
		//TODO: exception
		return 0;
	}*/

	double norm = 0;
	//	if(nParts==2){
	//norm by numerical integral
	//		norm = pPIF_->integral(minPar);
	//	}else if(nParts==3){
	//norm by phasespace monte-carlo
	if(pPHSP_){
		for(unsigned int phsp=0; phsp<nPhsp_; phsp++){
			Event theEvent(pPHSP_->getEvent(phsp));
			if(theEvent.getNParticles()!=3) continue;
			dataPoint point(theEvent);
			double intens = 0;
			if(pPIF_){
				ParameterList intensL = pPIF_->intensity(point, minPar);
				intens = intensL.GetDoubleParameter(0)->GetValue();
			}else{
				//TODO: Exception
				intens=0;
			}
			if(intens>0) norm+=intens;
		}
		//norm/=nPHSPEvts;
		//norm*=pPIF_->volume()/2.;
		//norm=nEvents*log(norm);
		//savedNorm=norm;
	}else if(pPhspTree_){
      pPhspTree_->recalculate();
      std::shared_ptr<DoubleParameter> intensL = std::dynamic_pointer_cast<DoubleParameter>(pPhspTree_->head()->getValue());
      norm = intensL->GetValue();
    }else{
	  norm=pPIF_->integral(minPar);
	}
	//	}

	/*std::cout << std::endl << "ControlPar LH: " << std::endl;
  std::cout << "Events: " << nEvents << std::endl;
  for(unsigned int i=0; i<minPar.GetNDouble(); i++){
    std::cout << minPar.GetParameterValue(i) << " ";
  }
  std::cout << std::endl;*/
	//std::cout << "BLA1" << std::endl;

	double lh=0; //calculate LH:
	//	switch(nParts){ //TODO: other cases, better "x" description (selection of particles?)
	//	case 2:{
	//		for(unsigned int evt = 0; evt < nEvents; evt++){
	//			Event theEvent(pDIF_->getEvent(evt));

	// TODO: try read exceptions

	//			std::vector<double> x;
	//			std::vector<double> xSq;
	//			if( nParts != theEvent.getNParticles()) continue; //TODO: real event count?
	//
	//			const Particle &a(theEvent.getParticle(0));
	//			const Particle &b(theEvent.getParticle(1));
	//
	//			double masssq = 0;
	//			masssq += (pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2));
	//			xSq.push_back(sqrt(masssq));

	//double intens = pPIF_->intensity(x, minPar);
	//ParameterList intensL = pPIF_->intensity(x, minPar);
	//double intens = intensL.GetDoubleParameter(0).GetValue();
	//			double intens = 0;
	//			if(pPIF_){
	//				ParameterList intensL = pPIF_->intensity(xSq, minPar);
	//				intens = intensL.GetDoubleParameter(0)->GetValue();
	//			}else if(pFcnTree_){
	//				//actualize inv masses
	//				minPar.GetDoubleParameter("ma")->SetValue(x[0]);
	//				minPar.GetDoubleParameter("mb")->SetValue(x[1]);
	//				minPar.GetDoubleParameter("mc")->SetValue(x[2]);
	//				//calculate intensity
	//				pFcnTree_->recalculate();
	//				std::shared_ptr<DoubleParameter> intensL = std::dynamic_pointer_cast<DoubleParameter>(pFcnTree_->head()->getValue());
	//				intens = intensL->GetValue();
	//			}else{
	//				//TODO: Exception
	//				intens=0;
	//			}
	//			if(intens>0){
	//				lh -= (log(intens/norm));
	//			}
	//		}

	//		break;
	//	}
	//	case 3:{
	if(pDIF_ && pPIF_){
	  for(unsigned int evt = 0; evt < nEvts_; evt++){
		Event theEvent(pDIF_->getEvent(evt));
		dataPoint point(theEvent);

		double intens = 0;
		if(pPIF_){
			ParameterList intensL = pPIF_->intensity(point, minPar);
			intens = intensL.GetDoubleParameter(0)->GetValue();
		}else{
			//TODO: Exception
			intens=0;
		}
		if(intens>0){
			lh += log(intens);
			//				std::cout<<"m23sq="<<x[0]<< " m13sq="<<x[1]<< " intens="<<intens<< " lh="<<lh<<std::endl;
		}
	  }
	}else if(pEvtTree_){
      pEvtTree_->recalculate();
      std::shared_ptr<DoubleParameter> intensL = std::dynamic_pointer_cast<DoubleParameter>(pEvtTree_->head()->getValue());
      lh = intensL->GetValue();
	}
	//lh = nEvents/2.*(norm/(nPHSPEvts-1))*(norm/(nPHSPEvts-1)) - lh + nEvents*log10(norm/nPHSPEvts);
//	std::cout.precision(15);
//	std::cout<<"event LH="<<lh<<" "<<nEvents<< " "<<norm/nPHSPEvts<<std::endl;
//	std::cout<<"phase space volume: "<<phspVolume<<std::endl;
	lh = nEvts_*log(norm/nPhsp_*phspVolume) - lh ;
//	std::cout<<"LH="<<lh<<std::endl;
	//lh -= norm;
//	break;
	//	}
	//	default:{
	//TODO: exception "i dont know how to handle this data"
	//		break;
	//	}
	//	}//end switch-case

	//std::cout << "BLA2" << std::endl;

	//std::cout << "ControlPar list " << minPar.GetNDouble() <<std::endl;

//	if(lh>0) BOOST_LOG_TRIVIAL(error) << "MinLogLH: positive -log(L)!";
	return lh;
}
