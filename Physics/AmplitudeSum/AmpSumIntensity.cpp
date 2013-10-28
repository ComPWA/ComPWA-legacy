
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

/// Default Constructor (0x0)
AmpSumIntensity::AmpSumIntensity(const DPKinematics kin, AmplitudeSetup ini) :
_kin(kin),
totAmp("relBWsumAmplitude", "totAmp"),
ampSetup(ini)
{
	init();
}

AmpSumIntensity::AmpSumIntensity(const double inM, const double inBr, const double in1
		,const double in2, const double in3, AmplitudeSetup ini) :
								  _kin(inM, inBr, in1, in2, in3,"","",""),
								  totAmp("relBWsumAmplitude", "totAmp"),
								  ampSetup(ini)
{

	init();
}

void AmpSumIntensity::init(){

	for(std::vector<Resonance>::iterator reso=ampSetup.getResonances().begin(); reso!=ampSetup.getResonances().end(); reso++){
		Resonance tmp = (*reso);
		//setup RooVars
		namer.push_back(tmp.m_name);
		mr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("mass_"+tmp.m_name,tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) );
		gr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("width_"+tmp.m_name,tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) );
		rr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("rr_"+tmp.m_name,tmp.m_strength) ) );
		phir.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("phir_"+tmp.m_name,tmp.m_phase) ) );

		//		qr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("par1_"+tmp.m_nametmp.m_breakup_mom) ) );
		//		aj.push_back( std::shared_ptr<IntegerParameter> (new IntegerParameter(tmp.m_spin) ) );
		//		am.push_back( std::shared_ptr<IntegerParameter> (new IntegerParameter(tmp.m_m) ) );
		//		an.push_back( std::shared_ptr<IntegerParameter> (new IntegerParameter(tmp.m_n) ) );

		//		par1.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("par1_"+tmp.m_name,tmp.m_par1) ) );
		//		par2.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("par1_"+tmp.name,tmp.m_par2) ) );

		unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;

		//setup Dynamics
		unsigned int last = mr.size()-1;
		if(tmp.m_type=="relBW"){
			std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
					*mr[last], *gr[last], tmp.m_breakup_mom, subSys, tmp.m_spin,tmp.m_m,tmp.m_n) );
			totAmp.addBW(tmpbw, rr.at(last), phir.at(last));
		}
		else if(tmp.m_type=="flatte"){
			std::shared_ptr<AmpFlatteRes> tmpbw(new AmpFlatteRes(tmp.m_name.c_str(),
					*mr[last], *gr[last], tmp.m_breakup_mom, *par1[last], *par2[last],0.547853, 0.1396, subSys, tmp.m_spin,tmp.m_m,tmp.m_n) );
			//			tmpbw->setBarrierMass(0.547853,0.1396);//a_0->eta pi hidden channel
			totAmp.addBW(tmpbw, rr.at(last), phir.at(last));
		}
		else {
			std::cout<<"AmpSumIntensity: wrong type! Type specified: "<<tmp.m_type<<std::endl;
			exit(1);
		}
	}// end loop over resonances
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
	if(x.size()!=2) {
		std::cout<<"AmpSumIntensity: wrong size of phase space variables!"<<std::endl;
		exit(1);
	}
	dataPoint::instance()->setM(2,3,x[0]);
	dataPoint::instance()->setM(1,3,x[1]);
	return intensity(par);
}
const ParameterList AmpSumIntensity::intensity( ParameterList& par){
	//parameters varied by Minimization algorithm
	for(unsigned int i=0; i<nAmps; i++){
		rr[i]->SetValue(par.GetDoubleParameter(i)->GetValue());//free
		phir[i]->SetValue(par.GetDoubleParameter(nAmps+i)->GetValue());//fixed
	}

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

void AmpSumIntensity::printAmps(){
	cout<<"== Printing amplitudes with current(!) set of parameters:"<<endl;
	for(unsigned int i=0;i<nAmps;i++)
		cout<<namer[i]<<":	Amplitude: "<<rr[i]->GetValue()<<"+-"<<rr[i]->GetError()<<"	Phase: "<<phir[i]->GetValue()<<"+-"<<phir[i]->GetError()<<endl;
}
