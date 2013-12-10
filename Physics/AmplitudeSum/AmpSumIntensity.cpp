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

#include "Physics/Amplitude.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
//#include "Estimator/AmpFcn.cpp"
//#include "Optimizer/Minuit2/MinuitIF.hpp"
//#include "Estimator/MinLogLH/MinLogLH.hpp"
#include <stdlib.h>
#include "gsl/gsl_monte.h"
#include "gsl/gsl_monte_plain.h"
#include "gsl/gsl_monte_miser.h"
#include "gsl/gsl_monte_vegas.h"

#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"

#include "Core/PhysConst.hpp"

#include "boost/function.hpp"


AmpSumIntensity::AmpSumIntensity(const AmpSumIntensity& other) : nAmps(other.nAmps), _dpArea(other._dpArea),
_entries(other._entries), _normStyle(other._normStyle),maxVal(other.maxVal),ampSetup(other.ampSetup),totAmp(other.totAmp){
}

/// Default Constructor (0x0)
AmpSumIntensity::AmpSumIntensity(AmplitudeSetup ini, unsigned int entries, normalizationStyle ns) :
		//_kin(kin),
		totAmp("relBWsumAmplitude", "totAmp"),
		ampSetup(ini),
		_entries(entries),
		_normStyle(ns)
{
	init();
}

AmpSumIntensity::AmpSumIntensity(const double inM, const double inBr, const double in1,const double in2, const double in3,
		std::string nameM, std::string name1,std::string name2,std::string name3,
		AmplitudeSetup ini, unsigned int entries, normalizationStyle ns) :
		//										 _kin(inM, inBr, in1, in2, in3,nameM,name1,name2,name3),
												 totAmp("relBWsumAmplitude", "totAmp"),
												 ampSetup(ini),
												 _entries(entries),
												 _normStyle(ns)
{
	init();
}

void AmpSumIntensity::init(){

	_dpArea = DalitzKinematics::instance()->getDParea();
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
		else if(_normStyle==one) norm = sqrt(tmpbw->integral());
		else if(_normStyle==entries) norm = sqrt(_dpArea*tmpbw->integral()/_entries);

		tmpbw->SetNormalization(1/norm);
		std::cout<<"AmpSumIntensity: INFO: Normalization constant for "<<tmp.m_name<<": "<<1/norm<<std::endl;
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
		else if(_normStyle==one) norm = sqrt(tmpbw->integral());
		else if(_normStyle==entries) norm = sqrt(_dpArea*tmpbw->integral()/_entries);
		tmpbw->SetNormalization(1/norm);
		std::cout<<"AmpSumIntensity: INFO: Normalization constant for "<<tmp.m_name<<": "<<1/norm<<std::endl;
	}// end loop over resonancesFlatte

	nAmps=rr.size();
	integral();
	std::cout << "completed setup" << std::endl;
}

double AmpSumIntensity::evaluate(double x[], size_t dim) {
	if(dim!=2) return 0;
	//set data point: we assume that x[0]=m13 and x[1]=m23
	dataPoint2 point; point.setVal("m13sq",x[0]); point.setVal("m23sq",x[1]);
	double m12sq = DalitzKinematics::instance()->getThirdVariableSq(x[0],x[1]);
	if( !DalitzKinematics::instance()->isWithinDP(x[1],x[0],m12sq) ) return 0;//only integrate over phase space
	ParameterList res = intensity(point);
	double intens = *res.GetDoubleParameter(0);
	return intens;
}
double evalWrapperAmpSumIntensity(double* x, size_t dim, void* param) {
	/* We need a wrapper here because intensity() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	return static_cast<AmpSumIntensity*>(param)->evaluate(x,dim);
};

const double AmpSumIntensity::integral(ParameterList& par){
	setParameterList(par);
	return integral();
}
const double AmpSumIntensity::integral(){

	/*
	 * integration functionality was tested with a model with only one normalized amplitude.
	 * The integration result is equal to the amplitude coefficient^2.
	 */

	size_t dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	DalitzKinematics* kin = DalitzKinematics::instance();
	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	//	double xLimit_low[2] = {0,0};
	//	double xLimit_high[2] = {10,10};
	size_t calls = 100000;
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function G = {&evalWrapperAmpSumIntensity,dim, const_cast<AmpSumIntensity*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&G, xLimit_low, xLimit_high, 2, calls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	std::cout<<"AmpSumIntensity: INFO: Integration result for amplitude sum: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res<<std::endl;

	return res;
}
const ParameterList AmpSumIntensity::intensity(std::vector<double> point, ParameterList& par){
	setParameterList(par);
	dataPoint2 dataP; dataP.setVal("m23sq",point[0]); dataP.setVal("m13sq",point[1]);
	return intensity(dataP);
}
const ParameterList AmpSumIntensity::intensity(dataPoint2& point, ParameterList& par){
	setParameterList(par);
	return intensity(point);
}
const ParameterList AmpSumIntensity::intensity(dataPoint2& point){
	//	std::cout<<dataPoint::instance()->getMsq(2,3)<<" "<<dataPoint::instance()->getMsq(1,3)<<" "<<dataPoint::instance()->getMsq(1,2)<<std::endl;
	double AMPpdf = totAmp.evaluate(point);
	if(AMPpdf!=AMPpdf){
		std::cout << "Error AmpSumIntensity: Intensity is not a number!!" << std::endl;
		exit(1);
	}
	ParameterList result;
	result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("AmpSumResult",AMPpdf)));
	return result;
}
void AmpSumIntensity::setParameterList(ParameterList& par){
	//parameters varied by Minimization algorithm
	for(unsigned int i=0; i<nAmps; i++){
		rr[i]->SetValue(par.GetDoubleParameter(2*i)->GetValue());//free
		phir[i]->SetValue(par.GetDoubleParameter(2*i+1)->GetValue());//fixed
	}
	return;
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
