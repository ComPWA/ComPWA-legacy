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

#include "Estimator/SliceFit/SliceFit.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/ParameterList.hpp"

namespace ComPWA {

using DataReader::Data;
using Physics::AmplitudeSum::AmpSumIntensity;
using Physics::DPKinematics::DalitzKinematics;
using Optimizer::ControlParameter;

namespace Estimator {
namespace SliceFit {

SliceFit::SliceFit(std::shared_ptr<AmpSumIntensity> inPIF, std::shared_ptr<Data> inDIF, ParameterList& inPar, unsigned int startEvent, unsigned int nEvents)
: pPIF_(inPIF), pDIF_(inDIF), nEvts_(0), nPhsp_(0), nStartEvt_(startEvent), par_(inPar), nUseEvt_(nEvents),nBins_(200){
	phspVolume = Kinematics::instance()->getPhspVolume();
	nEvts_ = pDIF_->getNEvents();
	if(startEvent+nUseEvt_<nEvts_) nUseEvt_ = nEvts_-startEvent;
	init();
}

SliceFit::SliceFit(std::shared_ptr<AmpSumIntensity> inPIF, std::shared_ptr<Data> inDIF, std::shared_ptr<Data> inPHSP, ParameterList& inPar, unsigned int startEvent, unsigned int nEvents)
: pPIF_(inPIF), pDIF_(inDIF), pPHSP_(inPHSP), nEvts_(0), nPhsp_(0), nStartEvt_(startEvent), par_(inPar), nUseEvt_(nEvents),nBins_(200){
	phspVolume = Kinematics::instance()->getPhspVolume();
	nEvts_ = pDIF_->getNEvents();
	nPhsp_ = inPHSP->getNEvents();
	if(!(startEvent+nUseEvt_<=nEvts_)) nUseEvt_ = nEvts_-startEvent;
	if(!(startEvent+nUseEvt_<=nPhsp_)) nUseEvt_ = nPhsp_-startEvent;
	init();
}

void SliceFit::init(){

	//aSlice_ = new TH1D("aSlice","Dummy",nBins_,0,10);
	//theAmp_ = new TH1D("theAmp","Dummy",nBins_,0,10);

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());

	M = 3.096916; // GeV/c² (J/psi+)
	Br = 0.000093; // GeV/c² (width)
	m1 = 0.; // GeV/c² (gamma)
	m2 = 0.139570; // GeV/c² (pi)
	m3 = 0.139570; // GeV/c² (pi)
	PI = 3.14159; // m/s

	Double_t m23_min = (kin->m2+kin->m3);
	Double_t m23_max = (kin->M-kin->m1);
	Double_t m13_min = (kin->m1+kin->m3);
	Double_t m13_max = (kin->M-kin->m2);
	Double_t m12_min = (kin->m1+kin->m2);
	Double_t m12_max = (kin->M-kin->m3);

	dalitzPlot_ = new TH2D("SliceDalitz","SliceDalitz", //TODO: variable binning
			nBins_,m23_min*m23_min,m23_max*m23_max, nBins_,m13_min*m13_min,m13_max*m13_max);

	//fill dalitz
	double m23sq, m13sq, m12sq;

	for(unsigned int i=0; i<pDIF_->getNEvents();i++){
		Event& event(pDIF_->getEvent(i));

		//if(!myReader.getEvent(i, event)) continue; TODO: try exception
		if(!event.getNParticles() == 3) continue;
		const Particle &a(event.getParticle(0));
		const Particle &b(event.getParticle(1));
		const Particle &c(event.getParticle(2));
		m12sq = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);
		m13sq = pow(a.E+c.E,2) - pow(a.px+c.px ,2) - pow(a.py+c.py ,2) - pow(a.pz+c.pz ,2);
		m23sq = pow(b.E+c.E,2) - pow(b.px+c.px ,2) - pow(b.py+c.py ,2) - pow(b.pz+c.pz ,2);

		dalitzPlot_->Fill(m23sq ,m13sq);
	}


}

std::shared_ptr<ControlParameter> SliceFit::createInstance(std::shared_ptr<AmpSumIntensity> inPIF, std::shared_ptr<Data> inDIF, ParameterList& inPar, unsigned int startEvent, unsigned int nEvents){
	if(!instance_)
		instance_ = std::shared_ptr<ControlParameter>(new SliceFit(inPIF, inDIF, inPar, startEvent, nEvents));

	return instance_;
}

std::shared_ptr<ControlParameter> SliceFit::createInstance(std::shared_ptr<AmpSumIntensity> inPIF, std::shared_ptr<Data> inDIF, std::shared_ptr<Data> inPHSP, ParameterList& inPar, unsigned int startEvent, unsigned int nEvents){
	if(!instance_)
		instance_ = std::shared_ptr<ControlParameter>(new SliceFit(inPIF, inDIF, inPHSP, inPar, startEvent, nEvents));

	return instance_;
}

SliceFit::~SliceFit(){
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

double SliceFit::controlParameter(ParameterList& minPar){
	unsigned int nEvents = pDIF_->getNEvents();
	unsigned int nParts = ((Event)pDIF_->getEvent(0)).getNParticles();

	//check if able to handle this many particles
	if(nParts!=3){
		//TODO: exception
		return 0;
	}

	//delete aSlice_;
	//delete theAmp_;

	std::string name1="Slice ", name2="Slice ", name3="Slice ";
	name1+=std::to_string((long long unsigned int) whichSlice_)+" data";
	name2+=std::to_string((long long unsigned int) whichSlice_)+" modelSlice";
	name3+=std::to_string((long long unsigned int) whichSlice_)+" modelClassic";
	aSlice_ = std::shared_ptr<TH1D>(new TH1D(name1.c_str(),name1.c_str(),nBins_,0,10));
	theAmpSl_ = std::shared_ptr<TH1D>(new TH1D(name2.c_str(),name2.c_str(),nBins_,0,10));
	theAmpCl_ = std::shared_ptr<TH1D>(new TH1D(name3.c_str(),name3.c_str(),nBins_,0,10));

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());

	Double_t m23_sq_min = (kin->m2+kin->m3); m23_sq_min*=m23_sq_min;
	Double_t m23_sq_max = (kin->M-kin->m1); m23_sq_max*=m23_sq_max;
	Double_t m13_sq_min = (kin->m1+kin->m3); m13_sq_min*=m13_sq_min;
	Double_t m13_sq_max = (kin->M-kin->m2); m13_sq_max*=m13_sq_max;
	Double_t m12_sq_min = (kin->m1+kin->m2); m12_sq_min*=m12_sq_min;
	Double_t m12_sq_max = (kin->M-kin->m3); m12_sq_max*=m12_sq_max;

	//------------fit slice number whichSlice_
	double lh=0;
	//for(int i=0+(nBins_/20.); i<nBins_-(nBins_/10.); i++){//slices
	double m23 = dalitzPlot_->GetBinCenter(whichSlice_);
	double locmin = m13_sq_min_constr(m23), locmax = m13_sq_max_constr(m23);


	for(unsigned int j=0; j<nBins_; j++){//slice values
		double m13 = dalitzPlot_->GetYaxis()->GetBinCenter(j);
		if(m13<locmin || m13>locmax) continue;

		double val = dalitzPlot_->GetBinContent(whichSlice_,j);

		dataPoint point;
		point.setVal("m23sq",m23); point.setVal("m13sq",m13);//point.push_back(m12sq);

		std::complex<double> reso[2];
		//reso[2]=std::complex<double>(minPar.GetParameterValue(5),0.);
		reso[0]=std::complex<double>(minPar.GetParameterValue(1),minPar.GetParameterValue(2));
		reso[1]=std::complex<double>(minPar.GetParameterValue(3),minPar.GetParameterValue(4));

		double amp = std::fabs(minPar.GetParameterValue(0))*pPIF_->sliceIntensity(point, par_,reso, 2);
		double ampCl = (double) (0.0725*(pPIF_->intensity(point, par_).GetParameterValue(0)));

		aSlice_->SetBinContent(j,val);
		theAmpSl_->SetBinContent(j,amp);
		theAmpCl_->SetBinContent(j,ampCl);

		double delta = (amp-(val*std::log(amp)));
		if(delta==delta)
			lh+=delta;//(val-amp)*(val-amp);///(double)nBins_; model_value - data_point->z * log(model_value);
		else{
			std::cout << "log( " << amp <<" ), point: " << m23 << ";"<<m13<<std::endl;
		}
	}//end slice-bins

	//}

	return lh;


	//TH1D *slice;
	//TF1 *slicefun;
	/*   TString th1name = "slice"; th1name+=i;
    slice = new TH1D(th1name,"m23 slice", nBins_,m13_sq_min,m13_sq_max);
    slice->SetTitle("");
    //slice->GetYaxis()->SetTitle("Magnitude");
    slice->GetXaxis()->SetTitle("m_{13}^{2} / GeV^{2}");
    double m23 = dalitzPlot_->GetBinCenter(i);
   // cout << endl << "Slice: " << i << "\t m23 = " << m23 << endl;
    double locmin = m13_sq_min_constr(m23), locmax = m13_sq_max_constr(m23);

    slicefun = new TF1("fitslice",&func_,locmin,locmax,7,"FitFuncObject");//TF1 * f = new TF1("f",fobj,0,1,npar,"MyFunctionObject");

    //slicefun->FixParameter(1,0.);//fix one parameter for first phase
    slicefun->FixParameter(0,m23);
    //slicefun->SetParameter(1,100);
    slicefun->SetParLimits(1,10,1000); slicefun->FixParameter(1,140);
    slicefun->SetParLimits(2,-10,10);
    slicefun->SetParLimits(3,-10,10);
    slicefun->SetParLimits(4,-10,10);
    slicefun->SetParLimits(5,-10,10);
    slicefun->SetParLimits(6,-10,10); slicefun->FixParameter(6,1.);

    for(unsigned int j=0; j<nBins_; j++){
      slice->SetBinContent(j,dalitzPlot_->GetBinContent(i,j));
    }//end slice-bins
	 */
	//if(i==32){
	//  reso1 = new TH1D(*slice);
	//  reso1->Fit(slicefun);
	//}else if(i==64){
	//  reso2 = new TH1D(*slice);
	//  reso2->Fit(slicefun);
	//slicefunAMP = new TF1("fitsliceAMP",fitsliceAMP,locmin,locmax,6);
	//slicefunAMP->FixParameter(0,m23);
	//slicefunAMP->SetParLimits(1,0,100000);
	//reso2->Fit(slicefunAMP);
	//}

	//slice->Fit(slicefun);
	/*
    slice->SetFillColor(kGray+2);
    slice->DrawCopy("hbar");
    canS->RedrawAxis();
    canS->Print("ampAllSliceFits.ps");
    if(gifCnt==4){
      gifCnt=0;
      canS->Print("ampAllSliceFits.gif+50");
    }
    gifCnt++;
   // if(meanpars[1]>1000) meanpars[1]=0;
    //if(meanpars[3]>1000) meanpars[3]=0;
   // if(meanpars[5]>1000) meanpars[5]=0;
    //if(meanpars[7]>1000) meanpars[5]=0;


    for(unsigned int j=1; j<3; j++){//(SpinsToCheck.size()+1); j++){
      double tmpx=0, tmpy=0, tmpxe=0, tmpye=0;
      tmpx=(slice->GetFunction("fitslice"))->GetParameter(2*j);
      tmpxe=(slice->GetFunction("fitslice"))->GetParError(2*j);
      tmpy=(slice->GetFunction("fitslice"))->GetParameter(2*j+1);
      tmpye=(slice->GetFunction("fitslice"))->GetParError(2*j+1);

      RooComplex eiph (tmpx, tmpy);

     // cout << endl << "x1 = " << tmpx << "\t y1 = " << tmpy << "\t phase = " << atan2(tmpx,tmpy) << endl << endl;

      double phi = atan2(tmpx,tmpy), phie=(tmpx*tmpye-tmpy*tmpxe)/eiph.abs()/eiph.abs(), abse=(tmpx*tmpxe+tmpy*tmpye)/eiph.abs();
      while(phi<0 && j==1){
        phi+=2*PI;
      };
     // while(phi>2*PI){
     //   phi-=2*PI;
     // };
      if(abse<0) abse*=-1.;
      if(phie<0) phie*=-1.;
      if(eiph.abs()<1.5 && abse){
        (par_r[j-1]).SetPoint(i-(NumBins/20.),m23,eiph.abs());
        (par_phi[j-1]).SetPoint(i-(NumBins/20.),m23,phi);
        (par_r[j-1]).SetPointError(i-(NumBins/20.),0,abse);
        (par_phi[j-1]).SetPointError(i-(NumBins/20.),0,phie);
      }
    }

    par_r0.SetPoint(i-(NumBins/20.),m23,(slice->GetFunction("fitslice"))->GetParameter(1));
    par_r0.SetPointError(i-(NumBins/20.),0,(slice->GetFunction("fitslice"))->GetParError(1));

    par_if.SetPoint(i-(NumBins/20.),m23,(slice->GetFunction("fitslice"))->GetParameter(6));
    par_if.SetPointError(i-(NumBins/20.),0,(slice->GetFunction("fitslice"))->GetParError(6));

    //cout << "x1 = " << (slice->GetFunction("fitslice"))->GetParameter(0) << "\t Err = " << (slice->GetFunction("fitslice"))->GetParError(0) << endl;
    //cout << "y1 = " << (slice->GetFunction("fitslice"))->GetParameter(1) << "\t Err = " << (slice->GetFunction("fitslice"))->GetParError(4) << endl;

    delete slice;
  }//end slices

  return 0;*/
}

} /* namespace SliceFit */
} /* namespace Estimator */
} /* namespace ComPWA */
