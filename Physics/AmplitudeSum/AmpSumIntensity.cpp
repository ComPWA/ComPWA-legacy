
#include <vector>
#include <memory>

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooCmdArg.h"
#include "RooMsgService.h"
#include "RooGlobalFunc.h"
#include "RooCFunction1Binding.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooComplex.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooRealProxy.h"

#include "Physics/Amplitude.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Estimator/AmpFcn.cpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"

#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"

/// Default Constructor (0x0)
AmpSumIntensity::AmpSumIntensity(const DPKinematics kin, AmplitudeSetup ini) :
_kin(kin),
dpPoint(_kin),
ma("ma", "mass", _kin.m23_min, _kin.m23_max),
mb("mb", "mass", _kin.m13_min, _kin.m13_max),
mc("mc", "mass", _kin.m12_min, _kin.m12_max),
totAmp("relBWsumAmplitude", "totAmp"),
ampSetup(ini)
{
	init();
}

AmpSumIntensity::AmpSumIntensity(const double inM, const double inBr, const double in1
		,const double in2, const double in3, AmplitudeSetup ini) :
				  _kin(inM, inBr, in1, in2, in3,"","",""),
				  dpPoint(_kin),
				  ma("ma", "mass", _kin.m23_min, _kin.m23_max),
				  mb("mb", "mass", _kin.m13_min, _kin.m13_max),
				  mc("mc", "mass", _kin.m12_min, _kin.m12_max),
				  totAmp("relBWsumAmplitude", "totAmp"),
				  ampSetup(ini)
{

	init();
}

void AmpSumIntensity::init(){
	//AmpSumOfAmplitudes totAmp23("rbwAmp23", "totAmp");

	//unsigned int numReso = ini.getResonances().size();
	for(std::vector<Resonance>::iterator reso=ampSetup.getResonances().begin(); reso!=ampSetup.getResonances().end(); reso++){
		Resonance tmp = (*reso);
		//setup RooVars
		namer.push_back(tmp.m_name);
		mr.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("m_{"+tmp.m_name+"}").c_str(), "resMass", tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) );
		qr.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("q_{"+tmp.m_name+"}").c_str(), "breakupMom", tmp.m_breakup_mom) ) );
		gr.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("g_{"+tmp.m_name+"}").c_str(), "resWidth", tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) );
		rr.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("r_{"+tmp.m_name+"}").c_str(), "amplitude", tmp.m_strength) ) );
		phir.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("phi_{"+tmp.m_name+"}").c_str(), "phase", tmp.m_phase) ) );

		aj.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("j_{"+tmp.m_name+"}").c_str(), "j", tmp.m_spin) ) );
		am.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("m_{"+tmp.m_name+"}").c_str(), "m", tmp.m_m) ) );
		an.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("n_{"+tmp.m_name+"}").c_str(), "n", tmp.m_n) ) );

		par1.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("par1_{"+tmp.m_name+"}").c_str(), "n", tmp.m_par1) ) );
		par2.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("par2_{"+tmp.m_name+"}").c_str(), "n", tmp.m_par2) ) );

		/*
		 * SORT ORDER: internally particles are sort as M -> (m1 m2) m3, with m1 m2 forming the resonance and m3 is the
		 * spectator hadron.
		 */
		int subSys = tmp.m_daugtherA + tmp.m_daugtherB;
		const double* mass_first;//masses of decay products
		const double* mass_third;
		const double* mass_second;
		RooAbsReal* mass12;//dalitz plot variables
		RooAbsReal* mass23;
		RooAbsReal* mass13;
		if((tmp.m_daugtherA==2 && tmp.m_daugtherB==3) ){
			mass_first = &_kin.m2;
			mass_second = &_kin.m3;
			mass_third = &_kin.m1;
			mass12 = &ma;
			mass23 = &mb;
			mass13 = &mc;
		}else if((tmp.m_daugtherA==3 && tmp.m_daugtherB==2) ){
			mass_first = &_kin.m3;
			mass_second = &_kin.m2;
			mass_third = &_kin.m1;
			mass12 = &ma;
			mass23 = &mb;
			mass13 = &mc;
		}else if( (tmp.m_daugtherA==1 && tmp.m_daugtherB==3) ){
			mass_first = &_kin.m1;
			mass_second = &_kin.m3;
			mass_third = &_kin.m2;
			mass12 = &mb;
			mass23 = &ma;
			mass13 = &mc;
		}else if( (tmp.m_daugtherA==3 && tmp.m_daugtherB==1) ){
			mass_first = &_kin.m3;
			mass_second = &_kin.m1;
			mass_third = &_kin.m2;
			mass12 = &mb;
			mass23 = &ma;
			mass13 = &mc;
		}else if( (tmp.m_daugtherA==1 && tmp.m_daugtherB==2) ){
			mass_first = &_kin.m1;//BUG SOMEWHERE for the case 1,2 or 2,1
			mass_second = &_kin.m2;
			mass_third = &_kin.m3;
			mass12 = &mc;
			mass23 = &ma;
			mass13 = &mb;
		}else if( (tmp.m_daugtherA==2 && tmp.m_daugtherB==1) ){
			mass_first = &_kin.m2;
			mass_second = &_kin.m1;
			mass_third = &_kin.m3;
			mass12 = &mc;
			mass23 = &ma;
			mass13 = &mb;
		}else{ //ignore resonance
			std::cout << "Resonance ignores due to ERROR!" << std::cout;
			mr.pop_back();
			qr.pop_back();
			gr.pop_back();
			rr.pop_back();
			phir.pop_back();
			aj.pop_back();
			am.pop_back();
			an.pop_back();
			continue;
		}

		//setup Dynamics and Angular Distribution
		unsigned int last = mr.size()-1;
		if(tmp.m_type=="relBW"){
//			std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
//					tmp.m_name.c_str(), *point, *mr[last], *gr[last], *qr[last], 1, tmp.m_spin,tmp.m_m,tmp.m_n) );
			std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
					tmp.m_name.c_str(),*mass23, *mass12, *mr[last], *gr[last], *qr[last], 1, tmp.m_spin,tmp.m_m,tmp.m_n) );
			tmpbw->setDecayMasses(*mass_first,*mass_second,*mass_third, _kin.M);
			rbw.push_back(tmpbw);
			//			std::shared_ptr<AmpWigner> tmpWigner(new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
			//					 *mass23,*mass12, tmp.m_spin, tmp.m_m, tmp.m_n)) ;
			//			tmpWigner->setDecayMasses(M, *mass_first, *mass_second , *mass_third);
			//			angd.push_back( tmpWigner );
			totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last));
		}
		else if(tmp.m_type=="flatte"){
			std::shared_ptr<AmpFlatteRes> tmpbw(new AmpFlatteRes(tmp.m_name.c_str(),
					tmp.m_name.c_str(),*mass23, *mass12, *mr[last], *gr[last], *qr[last], *par1[last], *par2[last], 1, tmp.m_spin,tmp.m_m,tmp.m_n) );
			tmpbw->setDecayMasses(*mass_first,*mass_second,*mass_third, _kin.M);
			tmpbw->setBarrierMass(0.547853,0.1396);//a_0->eta pi hidden channel
			rbw.push_back(tmpbw);
			//			std::shared_ptr<AmpWigner> tmpWigner(new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
			//					 *mass23,*mass12, tmp.m_spin, tmp.m_m, tmp.m_n)) ;
			//			tmpWigner->setDecayMasses(M, *mass_first, *mass_second , *mass_third);
			//			angd.push_back( tmpWigner );
			totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last));
		}
		else continue;
	}
	nAmps=rr.size();
	std::cout << "completed setup" << std::endl;
}
double AmpSumIntensity::getMaxVal(){

	std::cout<< "Calculating maximum value of function..."<<std::endl;
	std::shared_ptr<ControlParameter> esti( new AmpFcn(this));
	std::shared_ptr<Optimizer> opti(new MinuitIF(esti));
	ParameterList parList;
	//		parList.AddParameter(DoubleParameter((m23_max-m23_min)/2,m23_min,m23_max,0.001));//Start parameters
	//		parList.AddParameter(DoubleParameter((m13_max-m13_max)/2,m13_min,m13_max,0.001));//fit is sensitive to start values
	parList.AddParameter(DoubleParameter(1.1,_kin.m23_min,_kin.m23_max,0.001));//Start parameters
	parList.AddParameter(DoubleParameter(1.2,_kin.m13_min,_kin.m13_max,0.001));//fit is sensitive to start values
	//		std::vector<double> eeee;
	//		eeee.push_back(1.08);
	//		eeee.push_back(1.23);
	//		std::cout<<"intensity before fit: "<<this->intensity(eeee,eee)<<std::endl;
	ParameterList eee; this->fillStartParVec(eee);
	//  ParameterList startPar; this->fillStartParVec(startPar);
	//  for(unsigned int i=0;i<parList.GetNParameter();i++) cout<<parList.GetParameterValue(i)<<endl;
	//  for(unsigned int i=0;i<startPar.GetNParameter();i++) cout<<startPar.GetParameterValue(i)<<endl;
	maxVal = (-1)*opti->exec(parList);
	cout<<"maximum value: "<<maxVal<<endl;
	std::cout<<"maximum at: m23="<<ma.getVal()<<" m13="<<mb.getVal()<<" m12="<<mc.getVal()<<std::endl;
	return maxVal;
}

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

const double AmpSumIntensity::intensity(std::vector<double>& x, ParameterList& par){
	//TODO: check x exception
	if(x.size()!=2) return 0;

	ma.setVal(x[0]); mb.setVal(x[1]); //mc.setVal(x[2]);
	mc.setVal(_kin.M*_kin.M+_kin.m1*_kin.m1+_kin.m2*_kin.m2+_kin.m3*_kin.m3-x[0]*x[0]-x[1]*x[1]);
	//    if( par.GetNDouble()>rr.size() ){
	//        std::cout << "Error: Parameterlist doesn't match model!!" << std::endl; //TODO: exception
	//      return 0;
	//    }

	//parameters varied by Minimization algorithm
	for(unsigned int i=0; i<nAmps; i++){
		rr[i]->setVal(par.GetDoubleParameter(i).GetValue());//free
		phir[i]->setVal(par.GetDoubleParameter(nAmps+i).GetValue());//fixed
	}

	//	std::cout<<" -- "<<x[0]<<" "<<x[1]<<std::endl;
	double AMPpdf = totAmp.evaluate();

	if(AMPpdf!=AMPpdf){
		//        std::cout << "Error: Intensity is not a number!!" << std::endl; //TODO: exception
		return 0;
	}

	return AMPpdf;
}

const bool AmpSumIntensity::fillStartParVec(ParameterList& outPar){
	if(outPar.GetNParameter())
		return false; //already filled, TODO: exception?

	for(unsigned int i=0; i<rr.size();i++){//adding amplitudes
		if(rr[i]->hasError()){ //TODO: check bounds
			outPar.AddParameter(DoubleParameter(rr[i]->getVal(), rr[i]->getMin(), rr[i]->getMax(), rr[i]->getError()));
		}else{
			outPar.AddParameter(DoubleParameter(rr[i]->getVal(), 0.1));
		}
	}
	for(unsigned int i=0; i<phir.size();i++){//adding phases
		if(phir[i]->hasError()){ //TODO: check bounds
			outPar.AddParameter(DoubleParameter(phir[i]->getVal(), phir[i]->getMin(), phir[i]->getMax(), phir[i]->getError()));
		}else{
			outPar.AddParameter(DoubleParameter(phir[i]->getVal(), 0.1));
		}
	}
	return true;
}

void AmpSumIntensity::printAmps(){
	cout<<"== Printing amplitudes with current(!) set of parameters:"<<endl;
	for(unsigned int i=0;i<nAmps;i++)
		cout<<namer[i]<<":	Amplitude: "<<rr[i]->getVal()<<"+-"<<rr[i]->getError()<<"	Phase: "<<phir[i]->getVal()<<"+-"<<phir[i]->getError()<<endl;
}
