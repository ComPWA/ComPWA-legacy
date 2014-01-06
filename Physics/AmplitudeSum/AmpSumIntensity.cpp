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
#include <ctime>

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
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
//#include "Physics/AmplitudeSum/AmpSumTree.hpp"

#include "Core/PhysConst.hpp"

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include "boost/function.hpp"
#include <boost/log/trivial.hpp>
using namespace boost::log;


AmpSumIntensity::AmpSumIntensity(const AmpSumIntensity& other) : nAmps(other.nAmps), _dpArea(other._dpArea),
		_entries(other._entries), _normStyle(other._normStyle),maxVal(other.maxVal),ampSetup(other.ampSetup),
		totAmp(other.totAmp), _calcMaxFcnVal(other._calcMaxFcnVal), _maxFcnVal(other._maxFcnVal){
}

AmpSumIntensity::AmpSumIntensity(AmplitudeSetup ini, normalizationStyle ns, unsigned int entries, double dpArea) :
										totAmp("relBWsumAmplitude", "totAmp"), ampSetup(ini),
										_entries(entries), _normStyle(ns), _calcNorm(1), _dpArea(dpArea),
										_calcMaxFcnVal(0)
{
	init();
}

AmpSumIntensity::AmpSumIntensity(AmplitudeSetup ini, unsigned int entries, double dpArea) :
										totAmp("relBWsumAmplitude", "totAmp"), ampSetup(ini),
										_entries(entries), _normStyle(none), _calcNorm(0), _dpArea(dpArea),
										_calcMaxFcnVal(0)
{
	init();
}

AmpSumIntensity::AmpSumIntensity(const double inM, const double inBr, const double in1,const double in2, const double in3,
		std::string nameM, std::string name1,std::string name2,std::string name3,
		AmplitudeSetup ini, unsigned int entries, normalizationStyle ns, double dpArea) :
								totAmp("relBWsumAmplitude", "totAmp"), ampSetup(ini),
								_entries(entries), _normStyle(ns), _calcNorm(1),_dpArea(dpArea),
										_calcMaxFcnVal(0)
{
	init();
}

void AmpSumIntensity::init(){
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	if(_dpArea==-999) _dpArea = kin->getDParea();
	BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::init() number of Entries in dalitz plot set to: "<<_entries;
	BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::init() area of phase space: "<<_dpArea;

    //------------Setup Tree---------------------
    myTree = std::shared_ptr<FunctionTree>(new FunctionTree());

    //------------Setup Tree Pars---------------------
    std::shared_ptr<DoubleParameter> x = std::shared_ptr<DoubleParameter>(new DoubleParameter("x"));
    std::shared_ptr<DoubleParameter> m23 = std::shared_ptr<DoubleParameter>(new DoubleParameter("m23"));
    std::shared_ptr<DoubleParameter> m13 = std::shared_ptr<DoubleParameter>(new DoubleParameter("m13"));
    std::shared_ptr<DoubleParameter> m12 = std::shared_ptr<DoubleParameter>(new DoubleParameter("m12"));
    treePar = std::shared_ptr<ParameterList>(new ParameterList());
    treePar->AddParameter(x);
    treePar->AddParameter(m23);
    treePar->AddParameter(m13);
    treePar->AddParameter(m12);

    //----Strategies needed
    std::shared_ptr<MultAll> multStrat = std::shared_ptr<MultAll>(new MultAll());
    std::shared_ptr<AddAll> addStrat = std::shared_ptr<AddAll>(new AddAll());

    //----Add Top Node
    myTree->createHead("Amplitude", addStrat); //A=Sum{Resos}

    //----Add Resonances
	for(std::vector<Resonance>::iterator reso=ampSetup.getResonances().begin(); reso!=ampSetup.getResonances().end(); reso++){
		Resonance tmp = (*reso);
		if(!tmp.m_enable) continue;
		//setup RooVars
		namer.push_back(tmp.m_name);
		mr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("mass_"+tmp.m_name,tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) );
		gr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("width_"+tmp.m_name,tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) );
		rr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("rr_"+tmp.m_name,tmp.m_strength) ) );
		phir.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("phir_"+tmp.m_name,tmp.m_phase) ) );

		 //----Add Nodes
	    std::shared_ptr<BreitWignerStrategy> rbwStrat = std::shared_ptr<BreitWignerStrategy>(new BreitWignerStrategy(tmp.m_name));
	    std::shared_ptr<WignerDStrategy> angdStrat = std::shared_ptr<WignerDStrategy>(new WignerDStrategy(tmp.m_name));
	    unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;
		      unsigned int last = mr.size()-1;
		      myTree->createNode("Reso_"+tmp.m_name, multStrat, "Amplitude"); //Reso=BW*c*AD
		      myTree->createNode("RelBW_"+tmp.m_name, rbwStrat, "Reso_"+tmp.m_name); //BW
		      myTree->createLeaf("Intens_"+tmp.m_name, rr[last], "Reso_"+tmp.m_name); //c
		      myTree->createNode("AngD_"+tmp.m_name, angdStrat, "Reso_"+tmp.m_name); //AD
		      //BW Par
		      myTree->createLeaf("m0_"+tmp.m_name, mr[last]->GetValue(), "RelBW_"+tmp.m_name); //m0
		      myTree->createLeaf("x", x, "RelBW_"+tmp.m_name); //x
		      myTree->createLeaf("m23", m23, "RelBW_"+tmp.m_name); //ma TODO
		      myTree->createLeaf("m13", m13, "RelBW_"+tmp.m_name); //mb
		      myTree->createLeaf("spin_"+tmp.m_name, tmp.m_spin, "RelBW_"+tmp.m_name); //spin
		      myTree->createLeaf("d_"+tmp.m_name,  tmp.m_mesonRadius, "RelBW_"+tmp.m_name); //d
		      myTree->createLeaf("resWidth_"+tmp.m_name, gr[last]->GetValue(), "RelBW_"+tmp.m_name); //resWidth
		      //AD Par
		      myTree->createLeaf("m23", m23, "AngD_"+tmp.m_name); //ma TODO
		      myTree->createLeaf("m13", m13, "AngD_"+tmp.m_name); //mb
		      myTree->createLeaf("m12", m12, "AngD_"+tmp.m_name); //mc
		      myTree->createLeaf("M", kin->M, "AngD_"+tmp.m_name); //M
		      myTree->createLeaf("m1", kin->m1, "AngD_"+tmp.m_name); //m1
		      myTree->createLeaf("m2", kin->m2, "AngD_"+tmp.m_name); //m2
		      myTree->createLeaf("m3", kin->m3, "AngD_"+tmp.m_name); //m3
		     // unsigned int _subSysFlag = Double_t(paras.GetParameterValue("subSysFlag_"+name));
		      myTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "AngD_"+tmp.m_name); //subSysFlag
		      myTree->createLeaf("spin_"+tmp.m_name,tmp.m_spin, "AngD_"+tmp.m_name); //spin
		      myTree->createLeaf("m_"+tmp.m_name, tmp.m_m, "AngD_"+tmp.m_name); //OutSpin 1
		      myTree->createLeaf("n_"+tmp.m_name, tmp.m_n, "AngD_"+tmp.m_name); //OutSpin 2


		//setup Dynamics
		//unsigned int last = mr.size()-1;
		std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
				*mr[last], *gr[last], tmp.m_mesonRadius, subSys, tmp.m_spin,tmp.m_m,tmp.m_n) );
		totAmp.addBW(tmpbw, rr.at(last), phir.at(last));

		//setting normalization between amplitudes
		double norm=tmp.m_norm;
		if(norm<0 || _calcNorm) {//recalculate normalization
			norm = normReso(tmpbw);
			reso->m_norm = norm;//updating normalization
		}
		tmpbw->SetNormalization(1/norm);
	}// end loop over resonances


	for(std::vector<ResonanceFlatte>::iterator reso=ampSetup.getResonancesFlatte().begin(); reso!=ampSetup.getResonancesFlatte().end(); reso++){
		ResonanceFlatte tmp = (*reso);
		if(!tmp.m_enable) continue;
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

		double norm=tmp.m_norm;
		if(norm<0 || _calcNorm) {//recalculate normalization
			norm = normReso(tmpbw);
			reso->m_norm = norm;//updating normalization
		}
		tmpbw->SetNormalization(1/norm);
	}// end loop over resonancesFlatte

	ampSetup.save(ampSetup.getFilePath());//save updated information to input file
	nAmps=rr.size();
	if(_calcNorm) integral();
	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity: completed setup!";
}
double AmpSumIntensity::getMaxVal(){
	if(!_calcMaxFcnVal) calcMaxVal();
	return _maxFcnVal;
}
double AmpSumIntensity::getMaxVal(ParameterList& par){
	calcMaxVal(par);
	return _maxFcnVal;
}
void AmpSumIntensity::calcMaxVal(ParameterList& par){
	setParameterList(par);
	return calcMaxVal();
}
void AmpSumIntensity::calcMaxVal(){
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double maxVal=0;
	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	double maxM23=-999; double maxM13=-999;
	boost::minstd_rand rndGen2(123);
	boost::uniform_real<> uni_dist13(xLimit_low[0],xLimit_high[0]);
	boost::uniform_real<> uni_dist23(xLimit_low[1],xLimit_high[1]);
	boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uni13(rndGen2, uni_dist13);
	boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uni23(rndGen2, uni_dist23);
	for(unsigned int i=0; i<40000; i++){
		double m23sq=uni23(); double m13sq=uni13();
		dataPoint point; point.setVal("m13sq",m13sq); point.setVal("m23sq",m23sq);
		if( !kin->isWithinPhsp(point) ) { if(i>0) i--; continue; }//only integrate over phase space
		ParameterList res = intensity(point);
		double intens = *res.GetDoubleParameter(0);
		if(intens>maxVal) {
			maxM23=m23sq; maxM13=m13sq;
			maxVal=intens;
		}
	}
	_maxFcnVal=maxVal;
	_calcMaxFcnVal=1;
	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity::calcMaxVal() calculated maximum of amplitude: "
			<<_maxFcnVal<<" at m23sq="<<maxM23<<"/m13sq="<<maxM13;
	return ;
}
double AmpSumIntensity::normReso(std::shared_ptr<AmpAbsDynamicalFunction> amp){
	double norm;
	if(_normStyle==none) norm=1;
	else if(_normStyle==one) norm = sqrt(amp->integral());
	else if(_normStyle==entries) norm = sqrt(_dpArea*amp->integral()/_entries);
	BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::normRes Normalization constant for "
			<<amp->GetName()<<": "<<1/norm;
	return norm;

}
double AmpSumIntensity::evaluate(double x[], size_t dim) {
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//set data point: we assume that x[0]=m13 and x[1]=m23
	dataPoint point; point.setVal("m13sq",x[0]); point.setVal("m23sq",x[1]);
	double m12sq = kin->getThirdVariableSq(x[0],x[1]);
	if( !kin->isWithinPhsp(point) ) return 0;//only integrate over phase space
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
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
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
	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity::integrate() Integration result for amplitude sum: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}
const ParameterList AmpSumIntensity::intensity(std::vector<double> point, ParameterList& par){
	setParameterList(par);
	dataPoint dataP; dataP.setVal("m23sq",point[0]); dataP.setVal("m13sq",point[1]);
	return intensity(dataP);
}
const ParameterList AmpSumIntensity::intensity(dataPoint& point, ParameterList& par){
	setParameterList(par);
	return intensity(point);
}
const ParameterList AmpSumIntensity::intensity(dataPoint& point){
	//	std::cout<<dataPoint::instance()->getMsq(2,3)<<" "<<dataPoint::instance()->getMsq(1,3)<<" "<<dataPoint::instance()->getMsq(1,2)<<std::endl;
	double AMPpdf = totAmp.evaluate(point);
	if(AMPpdf!=AMPpdf){
		BOOST_LOG_TRIVIAL(error)<<"Error AmpSumIntensity: Intensity is not a number!!";
		exit(1);
	}
	ParameterList result;
	result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("AmpSumResult",AMPpdf)));
	return result;
}
std::shared_ptr<FunctionTree> AmpSumIntensity::functionTree(ParameterList& outPar) {
  if(outPar.GetNParameter()>0) return std::shared_ptr<FunctionTree>();
  fillStartParVec(outPar);
  outPar.AddParameter(treePar->GetDoubleParameter("x"));
  outPar.AddParameter(treePar->GetDoubleParameter("m23"));
  outPar.AddParameter(treePar->GetDoubleParameter("m13"));
  outPar.AddParameter(treePar->GetDoubleParameter("m12"));

  return myTree;
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
	BOOST_LOG_TRIVIAL(info)<<"== Printing amplitudes with current(!) set of parameters:";
	for(unsigned int i=0;i<nAmps;i++)
		BOOST_LOG_TRIVIAL(info)<<namer[i]<<":	Amplitude: "<<rr[i]->GetValue()<<"+-"<<rr[i]->GetError()<<"	Phase: "<<phir[i]->GetValue()<<"+-"<<phir[i]->GetError();
	return o.str();
}
