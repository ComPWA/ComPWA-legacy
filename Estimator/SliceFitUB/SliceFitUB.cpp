//-------------------------------------------------------------------------------
// Copyright (c) 2013 michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     michel - initial API and implementation
//-------------------------------------------------------------------------------
#include <sstream>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "Estimator/SliceFitUB/SliceFitUB.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/ParameterList.hpp"

namespace ComPWA {

using Physics::AmplitudeSum::AmpSumIntensity;
using Physics::DPKinematics::DalitzKinematics;

namespace Estimator {
namespace SliceFitUB {

SliceFitUB::SliceFitUB(std::shared_ptr<AmpSumIntensity> inPIF, std::shared_ptr<DataReader::Data> inDIF, ParameterList& inPar, unsigned int startEvent, unsigned int nEvents, unsigned int nBins, unsigned int nF0, unsigned int nF2)
: pPIF_(inPIF), pDIF_(inDIF), nEvts_(0), nPhsp_(0), nStartEvt_(startEvent), par_(inPar), nUseEvt_(nEvents),nBins_(nBins),nF0_(nF0),nF2_(nF2){
//phspVolume = Kinematics::instance()->getPhspVolume();
nEvts_ = pDIF_->getNEvents();
if(startEvent+nUseEvt_<nEvts_) nUseEvt_ = nEvts_-startEvent;
init();
}

SliceFitUB::SliceFitUB(std::shared_ptr<AmpSumIntensity> inPIF, std::shared_ptr<DataReader::Data> inDIF, std::shared_ptr<DataReader::Data> inPHSP, ParameterList& inPar, unsigned int startEvent, unsigned int nEvents, unsigned int nBins, unsigned int nF0, unsigned int nF2)
: pPIF_(inPIF), pDIF_(inDIF), pPHSP_(inPHSP), nEvts_(0), nPhsp_(0), nStartEvt_(startEvent), par_(inPar), nUseEvt_(nEvents),nBins_(nBins),nF0_(nF0),nF2_(nF2){
//phspVolume = Kinematics::instance()->getPhspVolume();
nEvts_ = pDIF_->getNEvents();
nPhsp_ = inPHSP->getNEvents();
if(!(startEvent+nUseEvt_<=nEvts_)) nUseEvt_ = nEvts_-startEvent;
if(!(startEvent+nUseEvt_<=nPhsp_)) nUseEvt_ = nPhsp_-startEvent;
init();
}

void SliceFitUB::init(){

  //aSlice_ = new TH1D("aSlice","Dummy",nBins_,0,10);
  //theAmp_ = new TH1D("theAmp","Dummy",nBins_,0,10);

  ComPWA::Physics::DPKinematics::DalitzKinematics* kin =
		  dynamic_cast<ComPWA::Physics::DPKinematics::DalitzKinematics*>(
				  Kinematics::instance()
  );

  M = 3.096916; // GeV/c² (J/psi+)
  Br = 0.000093; // GeV/c² (width)
  m1 = 0.; // GeV/c² (gamma)
  m2 = 0.139570; // GeV/c² (pi)
  m3 = 0.139570; // GeV/c² (pi)
  PI = 3.14159; // m/s

  double m23_min = (kin->m2+kin->m3);
  double m23_max = (kin->GetMotherMass()-kin->m1);
  double m13_min = (kin->m1+kin->m3);
  double m13_max = (kin->GetMotherMass()-kin->m2);
  double m12_min = (kin->m1+kin->m2);
  double m12_max = (kin->GetMotherMass()-kin->m3);

  // nBins_,m23_min*m23_min,m23_max*m23_max
  //double phspComplete = kin->getPhspVolume();
  double sliceWidth = (m23_max*m23_max-m23_min*m23_min)/nBins_;
  BOOST_LOG_TRIVIAL(debug) << "Slice width: " << sliceWidth << "\t of " << nBins_ << " bins";
  for(unsigned int i=0; i<nBins_; i++){
    slicedEvents_.push_back(std::vector<unsigned int>());
    slicedPhspEvt_.push_back(std::vector<unsigned int>());
    sliceMass_.push_back(m23_min*m23_min+(i+0.5)*(sliceWidth));
    slicedEvtMass_.push_back(std::vector<double>());

    phspVolume_.push_back(kin->getPhspVolumePart(m23_min*m23_min+(i)*(sliceWidth), m23_min*m23_min+(i+1)*(sliceWidth)));
  }

  double m23sq, m13sq, m12sq;

  for(unsigned int i=nStartEvt_; i<nUseEvt_+nStartEvt_;i++){
      Event& event(pDIF_->getEvent(i));

      //if(!myReader.getEvent(i, event)) continue; TODO: try exception
      if(!event.getNParticles() == 3) continue;
      const ComPWA::Particle &a(event.getParticle(0));
      const ComPWA::Particle &b(event.getParticle(1));
      const ComPWA::Particle &c(event.getParticle(2));
     // m12sq = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);
      m13sq = pow(a.E+c.E,2) - pow(a.px+c.px ,2) - pow(a.py+c.py ,2) - pow(a.pz+c.pz ,2);
      m23sq = pow(b.E+c.E,2) - pow(b.px+c.px ,2) - pow(b.py+c.py ,2) - pow(b.pz+c.pz ,2);

      //----separate events per slice
      int sliceID = floor((m23sq-m23_min*m23_min)/sliceWidth);
      if(sliceID<0 || sliceID>=nBins_)
        std::cout << "ERROR slicing event " << i <<" , ID & point: " << sliceID << " , " << m23sq << ";"<<m13sq<<std::endl;
      slicedEvents_[sliceID].push_back(i);
      slicedEvtMass_[sliceID].push_back(m23sq);

      //dalitzPlot_->Fill(m23sq ,m13sq);

  }

  for(unsigned int i=nStartEvt_; i<nUseEvt_+nStartEvt_;i++){
      Event& event(pPHSP_->getEvent(i));

      //if(!myReader.getEvent(i, event)) continue; TODO: try exception
      if(!event.getNParticles() == 3) continue;
      const ComPWA::Particle &a(event.getParticle(0));
      const ComPWA::Particle &b(event.getParticle(1));
      const ComPWA::Particle &c(event.getParticle(2));
     // m12sq = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);
      m13sq = pow(a.E+c.E,2) - pow(a.px+c.px ,2) - pow(a.py+c.py ,2) - pow(a.pz+c.pz ,2);
      m23sq = pow(b.E+c.E,2) - pow(b.px+c.px ,2) - pow(b.py+c.py ,2) - pow(b.pz+c.pz ,2);

      //----separate events per slice
      int sliceID = floor((m23sq-m23_min*m23_min)/sliceWidth);
      if(sliceID<0 || sliceID>=nBins_)
        std::cout << "ERROR slicing event " << i <<" , ID & point: " << sliceID << " , " << m23sq << ";"<<m13sq<<std::endl;
      slicedPhspEvt_[sliceID].push_back(i);

      //dalitzPlot_->Fill(m23sq ,m13sq);

  }

  for(unsigned int i=0; i<slicedEvents_.size(); i++){
    BOOST_LOG_TRIVIAL(debug) << "Slice " << i << "\t nEvents: " << slicedEvents_[i].size();
    if(slicedEvents_[i].size()>3){
      BOOST_LOG_TRIVIAL(debug) << "Three masses: " << slicedEvtMass_[i].at(0) << " " << slicedEvtMass_[i].at(1) << " " << slicedEvtMass_[i].at(2);
      BOOST_LOG_TRIVIAL(debug) << "in Slice: " << sliceMass_[i]-(sliceWidth/2.) << " to " << sliceMass_[i]+(sliceWidth/2.);
      BOOST_LOG_TRIVIAL(debug) << " ";
    }
  }

}

std::shared_ptr<Optimizer::ControlParameter> SliceFitUB::createInstance(std::shared_ptr<AmpSumIntensity> inPIF, std::shared_ptr<DataReader::Data> inDIF, ParameterList& inPar, unsigned int startEvent, unsigned int nEvents, unsigned int nBins, unsigned int nF0, unsigned int nF2){
    if(!instance_)
        instance_ = std::shared_ptr<Optimizer::ControlParameter>(new SliceFitUB(inPIF, inDIF, inPar, startEvent, nEvents, nBins, nF0, nF2));

    return instance_;
}

std::shared_ptr<Optimizer::ControlParameter> SliceFitUB::createInstance(std::shared_ptr<AmpSumIntensity> inPIF, std::shared_ptr<DataReader::Data> inDIF, std::shared_ptr<DataReader::Data> inPHSP, ParameterList& inPar, unsigned int startEvent, unsigned int nEvents, unsigned int nBins, unsigned int nF0, unsigned int nF2){
    if(!instance_)
        instance_ = std::shared_ptr<Optimizer::ControlParameter>(new SliceFitUB(inPIF, inDIF, inPHSP, inPar, startEvent, nEvents, nBins, nF0, nF2));

    return instance_;
}

SliceFitUB::~SliceFitUB(){
 // delete aSlice_;
 // delete theAmp_;
}

/*double SliceFit::fitsliceAMP(Double_t *x, Double_t *par){

  //double m12sq=M*M-x[0]-par[0];
 // if(m12sq<0)
  //  m12sq=0.0001;

  std::vector<double> point;
  point.push_back(par[0]); point.push_back(x[0]); //point.push_back(m12sq);

  std::complex<double> reso[3];
  reso[2]=std::complex<double>(par[6],0.);
  reso[0]=std::complex<double>(par[2],par[3]);
  reso[1]=std::complex<double>(par[4],par[5]);
  std::shared_ptr<AmpSumIntensity> amp = SliceFit::instance();
  double result = par[1]*amp->sliceIntensity(point, par_,reso, 3);
  //double result = par[1]*totAmp23.evaluate();

  if(result!=result) return 0;
  return result;
  //return totAmp23.evaluate();
}*/

double SliceFitUB::controlParameter(ComPWA::ParameterList& minPar){
  unsigned int nEvents = pDIF_->getNEvents();
  unsigned int nParts = ((Event)pDIF_->getEvent(0)).getNParticles();

  //check if able to handle this many particles
  if(nParts!=3){
      //TODO: exception
      return 0;
  }

  DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());

  double m23_min = (kin->m2+kin->m3);
  double m23_max = (kin->GetMotherMass()-kin->m1);
  double m13_min = (kin->m1+kin->m3);
  double m13_max = (kin->GetMotherMass()-kin->m2);
  double m12_min = (kin->m1+kin->m2);
  double m12_max = (kin->GetMotherMass()-kin->m3);

  std::string name1="Slice ", name2="Slice ", name3="Slice ";
  name1+=std::to_string(whichSlice_)+" data";
  name2+=std::to_string(whichSlice_)+" modelSlice";
  name3+=std::to_string(whichSlice_)+" modelClassic";
  aSlice_ = std::shared_ptr<TH1D>(new TH1D(name1.c_str(),name1.c_str(),nBins_,m13_min*m13_min,m13_max*m13_max));
  theAmpSl_ = std::shared_ptr<TH1D>(new TH1D(name2.c_str(),name2.c_str(),nBins_,m13_min*m13_min,m13_max*m13_max));
  theAmpCl_ = std::shared_ptr<TH1D>(new TH1D(name3.c_str(),name3.c_str(),nBins_,m13_min*m13_min,m13_max*m13_max));

  std::vector<unsigned int>& phspVec = slicedPhspEvt_[whichSlice_];
  std::vector<unsigned int>& dataVec = slicedEvents_[whichSlice_];

  std::complex<double> reso[2];
  //reso[2]=std::complex<double>(minPar.GetParameterValue(5),0.);
  reso[0]=std::complex<double>(
		  minPar.GetDoubleParameterValue(1),
		  minPar.GetDoubleParameterValue(2)
		  );
  reso[1]=std::complex<double>(
		  minPar.GetDoubleParameterValue(3),
		  minPar.GetDoubleParameterValue(4)
  );

  double norm = 0;
  if(pPHSP_){
      for(unsigned int phsp=0; phsp<phspVec.size(); phsp++){
          Event theEvent(pPHSP_->getEvent(phspVec[phsp]));
          if(theEvent.getNParticles()!=3) continue;
          dataPoint point(theEvent);
          double intens = 0;
          if(pPIF_){
              intens = pPIF_->sliceIntensity(
            		  point, par_,reso, 2,
					  std::fabs(minPar.GetDoubleParameterValue(0)),
					  nF0_,nF2_
              );//pPIF_->intensity(point, minPar);
              //intens = intensL.GetDoubleParameter(0)->GetValue();
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
  }//else{
  //    norm=pPIF_->integral(minPar);
  //}

  //norm = norm * phspVolume_[whichSlice_]/slicedPhspEvt_[whichSlice_].size();
  if(norm==0) norm=1;

  double lh=0;

  if(pDIF_ && pPIF_){
    for(unsigned int evt = 0; evt<dataVec.size(); evt++){
      Event theEvent(pDIF_->getEvent(dataVec[evt]));
      dataPoint point(theEvent);

      aSlice_->Fill(point.getVal("m13sq"));

      double intens = 0;
      if(pPIF_){
          intens = pPIF_->sliceIntensity(point, par_,reso, 2,
        		  std::fabs(minPar.GetDoubleParameterValue(0)),
				  nF0_,nF2_
          );//pPIF_->intensity(point, minPar);
          //intens = intensL.GetDoubleParameter(0)->GetValue();
      }else{
          //TODO: Exception
          intens=0;
      }
      if(intens>0){
        //lh += std::log(intens);
        lh += std::log(intens);
          //              std::cout<<"m23sq="<<x[0]<< " m13sq="<<x[1]<< " intens="<<intens<< " lh="<<lh<<std::endl;
      }

          //if(!evt)
          //  BOOST_LOG_TRIVIAL(debug) << "First Evt LH: " << lh;
      }
  }

  double sliceWidth = (m23_max*m23_max-m23_min*m23_min)/nBins_;
  double m23sq = (whichSlice_+0.5)*sliceWidth+m23_min*m23_min;  //id = floor((m23sq-m23_min*m23_min)/sliceWidth);
  double locmin = m13_sq_min_constr(m23sq), locmax = m13_sq_max_constr(m23sq);
  for(unsigned int bin=0; bin < nBins_; bin++){
    //?????? Fill Events in Histogram,  ...

    double m13sq = ((bin+0.5)*(m13_max*m13_max-m13_min*m13_min)/nBins_)+m13_min*m13_min;
    if(m13sq<locmin || m13sq>locmax) continue;

    //double val = dalitzPlot_->GetBinContent(whichSlice_,j);

    dataPoint point;
    point.setVal("m23sq",m23sq); point.setVal("m13sq",m13sq);//point.push_back(m12sq);

    std::complex<double> reso[2];
    //reso[2]=std::complex<double>(minPar.GetParameterValue(5),0.);
    reso[0]=std::complex<double>(
    		minPar.GetDoubleParameterValue(1),
			minPar.GetDoubleParameterValue(2)
    );
    reso[1]=std::complex<double>(
    		minPar.GetDoubleParameterValue(3),
			minPar.GetDoubleParameterValue(4)
    );

    double amp = pPIF_->sliceIntensity(
    		point, par_,reso, 2,
			std::fabs(minPar.GetDoubleParameterValue(0)),
			nF0_,nF2_
    );
    double ampCl = (double) (0.0725*(pPIF_->intensity(point).GetDoubleParameterValue(0)));

    theAmpSl_->SetBinContent(bin,amp);
    theAmpCl_->SetBinContent(bin,ampCl);
  }

  BOOST_LOG_TRIVIAL(debug) << "Data Term: " << lh << "\t Phsp Term (wo log): " << norm << " N " << slicedEvents_[whichSlice_].size() << " Np " << slicedPhspEvt_[whichSlice_].size() << " V " << phspVolume_[whichSlice_];
  lh = slicedEvents_[whichSlice_].size()*std::log(norm/slicedPhspEvt_[whichSlice_].size()*phspVolume_[whichSlice_]) - lh ;
  //lh = (norm/slicedPhspEvt_[whichSlice_].size()*phspVolume_[whichSlice_]) - lh ; //stefans note ext lh
  //lh = (-1)*((double)slicedEvents_[whichSlice_].size())/slicedEvents_[whichSlice_].size()*lh;
  //funcVal = sum(logLikelihoodAcc) + nmbEvt * sum(normFactorAcc); rootpwa
  return lh;
}

} /* namespace SliceFitUB */
} /* namespace Estimator */
} /* namespace ComPWA */
