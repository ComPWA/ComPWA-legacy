//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding flatte type resonance, removing root dependence
//-------------------------------------------------------------------------------

#include <vector>
#include <memory>

//#include "RooRealVar.h"
//#include "RooFormulaVar.h"
//#include "RooPlot.h"
//#include "RooCmdArg.h"
//#include "RooMsgService.h"
//#include "RooGlobalFunc.h"
//#include "RooCFunction1Binding.h"
//#include "RooGaussian.h"
//#include "RooAddPdf.h"
//#include "RooDataSet.h"
//#include "RooDataHist.h"
//#include "RooFitResult.h"
//#include "RooComplex.h"
//#include "RooAbsReal.h"
//#include "RooAbsArg.h"
//#include "RooRealProxy.h"

#include "Physics/Amplitude.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
//#include "Estimator/AmpFcn.cpp"
//#include "Optimizer/Minuit2/MinuitIF.hpp"
//#include "Estimator/MinLogLH/MinLogLH.hpp"

#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"

#include "Physics/DPKinematics/PhysConst.hpp"

#include "boost/function.hpp"
/// Default Constructor (0x0)
AmpSumIntensity::AmpSumIntensity(const DPKinematics kin, AmplitudeSetup ini, unsigned int entries, normalizationStyle ns) :
_kin(kin),
totAmp("relBWsumAmplitude", "totAmp"),
ampSetup(ini),
_entries(entries),
_normStyle(ns)
{
	init();
}

AmpSumIntensity::AmpSumIntensity(const double inM, const double inBr, const double in1
		,const double in2, const double in3, AmplitudeSetup ini, unsigned int entries, normalizationStyle ns) :
			 _kin(inM, inBr, in1, in2, in3,"","",""),
			 totAmp("relBWsumAmplitude", "totAmp"),
			 ampSetup(ini),
			 _entries(entries),
_normStyle(ns)
{
	init();
}

void AmpSumIntensity::init(){

	_dpArea = dataPoint::instance()->DPKin.getDParea();
	std::cout<<"AmpSumIntensity: INFO: number of Entries in dalitz plot set to: "<<_entries<<std::endl;
	std::cout<<"AmpSumIntensity: INFO: area of phase space: "<<_dpArea<<std::endl;

	for(std::vector<Resonance>::iterator reso=ampSetup.getResonances().begin(); reso!=ampSetup.getResonances().end(); reso++){
		Resonance tmp = (*reso);
		//setup RooVars
		namer.push_back(tmp.m_name);
		mr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("mass_"+tmp.m_name,tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) );
		gr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("width_"+tmp.m_name,tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) );
		rr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("rr_"+tmp.m_name,tmp.m_strength) ) );
		phir.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("phir_"+tmp.m_name,tmp.m_phase) ) );

		unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;

		//setup Dynamics
		unsigned int last = mr.size()-1;
		std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
				*mr[last], *gr[last], tmp.m_mesonRadius, subSys, tmp.m_spin,tmp.m_m,tmp.m_n) );
		totAmp.addBW(tmpbw, rr.at(last), phir.at(last));

		//setting normalization between amplitudes
		double norm=1;
		if(_normStyle==none) norm=1;
		else if(_normStyle==one) norm = tmpbw->integralNorm();
		else if(_normStyle==entries) norm = sqrt(_dpArea*tmpbw->integralNorm()/_entries);

		tmpbw->SetNormalization(1/norm);
		std::cout<<"AmpSumIntensity: INFO: Normalization constant for "<<tmp.m_name<<": "<<norm<<std::endl;
	}// end loop over resonances


	for(std::vector<ResonanceFlatte>::iterator reso=ampSetup.getResonancesFlatte().begin(); reso!=ampSetup.getResonancesFlatte().end(); reso++){
		ResonanceFlatte tmp = (*reso);
		//setup RooVars
		namer.push_back(tmp.m_name);
		mr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("mass_"+tmp.m_name,tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) );
		gr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("width_"+tmp.m_name,tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) );
		rr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("rr_"+tmp.m_name,tmp.m_strength) ) );
		phir.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("phir_"+tmp.m_name,tmp.m_phase) ) );
		DoubleParameter param1("coupling1_"+tmp.m_name,tmp.m_coupling);
		DoubleParameter param2("coupling2_"+tmp.m_name,tmp.m_couplingHidden);

		unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;

		//setup Dynamics
		unsigned int last = mr.size()-1;
		std::shared_ptr<AmpFlatteRes> tmpbw(new AmpFlatteRes(tmp.m_name.c_str(),
				*mr[last], *gr[last], tmp.m_mesonRadius, param1, param2, \
				PhysConst::instance()->getMass(tmp.m_hiddenParticle1),\
				PhysConst::instance()->getMass(tmp.m_hiddenParticle2),\
				subSys, tmp.m_spin,tmp.m_m,tmp.m_n) );
		totAmp.addBW(tmpbw, rr.at(last), phir.at(last));

		double norm=1;
		if(_normStyle==none) norm=1;
		else if(_normStyle==one) norm = tmpbw->integralNorm();
		else if(_normStyle==entries) norm = sqrt(_dpArea*tmpbw->integralNorm()/_entries);
		tmpbw->SetNormalization(1/norm);
		std::cout<<"AmpSumIntensity: INFO: Normalization constant for "<<tmp.m_name<<": "<<norm<<std::endl;
	}// end loop over resonancesFlatte

	nAmps=rr.size();
	std::cout << "completed setup" << std::endl;
}
//double AmpSumIntensity::getMaxVal(){

//	std::cout<< "Calculating maximum value of function..."<<std::endl;
//	std::shared_ptr<ControlParameter> esti( new AmpFcn(this));
//	std::shared_ptr<Optimizer> opti(new MinuitIF(esti));
//	ParameterList parList;
//	//		parList.AddParameter(DoubleParameter((m23_max-m23_min)/2,m23_min,m23_max,0.001));//Start parameters
//	//		parList.AddParameter(DoubleParameter((m13_max-m13_max)/2,m13_min,m13_max,0.001));//fit is sensitive to start values
//	parList.AddParameter(DoubleParameter(1.1,_kin.m23_min,_kin.m23_max,0.001));//Start parameters
//	parList.AddParameter(DoubleParameter(1.2,_kin.m13_min,_kin.m13_max,0.001));//fit is sensitive to start values
//	//		std::vector<double> eeee;
//	//		eeee.push_back(1.08);
//	//		eeee.push_back(1.23);
//	//		std::cout<<"intensity before fit: "<<this->intensity(eeee,eee)<<std::endl;
//	ParameterList eee; this->fillStartParVec(eee);
//	//  ParameterList startPar; this->fillStartParVec(startPar);
//	//  for(unsigned int i=0;i<parList.GetNParameter();i++) cout<<parList.GetParameterValue(i)<<endl;
//	//  for(unsigned int i=0;i<startPar.GetNParameter();i++) cout<<startPar.GetParameterValue(i)<<endl;
//	maxVal = (-1)*opti->exec(parList);
////	cout<<"maximum value: "<<maxVal<<endl;
////	std::cout<<"maximum at: m23="<<ma.getVal()<<" m13="<<mb.getVal()<<" m12="<<mc.getVal()<<std::endl;
//	return maxVal;
//	return 1;
//}

const double AmpSumIntensity::integral(ParameterList& par){
	/*double integral = 1;
    unsigned int nSteps = 1000000;
    double stepa = (ma.getMax()-ma.getMin())/(double)nSteps;
    double stepb = (mb.getMax()-mb.getMin())/(double)nSteps;

    //TODO: better approximation

    for(unsigned int k=1; k<nSteps; k++){
      integral += step*intensity((ma.getMin()+k*step), par);
    }*/

	return 1;
	//return totAmp.getNorm();//integral;
}
const ParameterList AmpSumIntensity::intensity(std::vector<double>& x, ParameterList& par){
	if(x.size()!=3) {
		std::cout<<"AmpSumIntensity: wrong size of phase space variables!"<<std::endl;
		exit(1);
	}
	dataPoint::instance()->setM(2,3,x[0]);
	dataPoint::instance()->setM(1,3,x[1]);
	dataPoint::instance()->setM(1,2,x[2]);
	return intensity(par);
}
const ParameterList AmpSumIntensity::intensity( ParameterList& par){
	//parameters varied by Minimization algorithm
	for(unsigned int i=0; i<nAmps; i++){
		rr[i]->SetValue(par.GetDoubleParameter(2*i)->GetValue());//free
		phir[i]->SetValue(par.GetDoubleParameter(2*i+1)->GetValue());//fixed
	}

	//	std::cout<<dataPoint::instance()->getMsq(2,3)<<" "<<dataPoint::instance()->getMsq(1,3)<<" "<<dataPoint::instance()->getMsq(1,2)<<std::endl;
	double AMPpdf = totAmp.evaluate();

	if(AMPpdf!=AMPpdf){
		std::cout << "Error AmpSumIntensity: Intensity is not a number!!" << std::endl;
		exit(1);
	}

	ParameterList result;
	result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("AmpSumResult",AMPpdf)));
	return result;
}

const bool AmpSumIntensity::fillStartParVec(ParameterList& outPar){
	if(outPar.GetNParameter())
		return false; //already filled, TODO: exception?
	for(unsigned int i=0; i<rr.size();i++){
		//add strength and phases of the used amplitudes
		outPar.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter(*rr[i])));
		outPar.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter(*phir[i])));
	}
	return true;
}

std::string AmpSumIntensity::printAmps(){
	std::stringstream o;
	o<<"== Printing amplitudes with current(!) set of parameters:"<<endl;
	for(unsigned int i=0;i<nAmps;i++)
		o<<namer[i]<<":	Amplitude: "<<rr[i]->GetValue()<<"+-"<<rr[i]->GetError()<<"	Phase: "<<phir[i]->GetValue()<<"+-"<<phir[i]->GetError()<<endl;
    return o.str();
}
