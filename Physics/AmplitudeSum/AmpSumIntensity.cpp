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
#include "Core/Functions.hpp"
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

AmpSumIntensity::AmpSumIntensity(AmplitudeSetup ini, normStyle ns, std::shared_ptr<Efficiency> eff, unsigned int entries, double dpArea) :
												totAmp("relBWsumAmplitude", "totAmp"), ampSetup(ini),
												_entries(entries), _normStyle(ns), _calcNorm(1), _dpArea(dpArea),
												_calcMaxFcnVal(0),eff_(eff)
{
	init();
}

AmpSumIntensity::AmpSumIntensity(AmplitudeSetup ini, std::shared_ptr<Efficiency> eff, unsigned int entries, double dpArea) :
												totAmp("relBWsumAmplitude", "totAmp"), ampSetup(ini),
												_entries(entries), _normStyle(none), _calcNorm(0), _dpArea(dpArea),
												_calcMaxFcnVal(0),eff_(eff)
{
	init();
}

AmpSumIntensity::AmpSumIntensity(const double inM, const double inBr, const double in1,const double in2, const double in3,
		std::string nameM, std::string name1,std::string name2,std::string name3,
		AmplitudeSetup ini, unsigned int entries, normStyle ns, double dpArea) :
												totAmp("relBWsumAmplitude", "totAmp"), ampSetup(ini),
												_entries(entries), _normStyle(ns), _calcNorm(1),_dpArea(dpArea),
												_calcMaxFcnVal(0),eff_(std::shared_ptr<Efficiency>(new UnitEfficiency()))
{
	init();
}

void AmpSumIntensity::init(){
	result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("AmpSumResult")));

	  DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	  _dpArea = kin->getPhspVolume();
	  //  std::cout<<kin->getPhspVolume()<<std::endl;
	  _calcNorm=1;
	  BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::init() number of Entries in dalitz plot set to: "<<_entries;
	  BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::init() area of phase space: "<<_dpArea;

	  for(std::vector<Resonance>::iterator reso=ampSetup.getResonances().begin(); reso!=ampSetup.getResonances().end(); reso++){
	       Resonance tmp = (*reso);
	       if(!tmp.m_enable) continue;
	       //setup RooVars
	       namer.push_back(tmp.m_name);
	       mr.insert(std::make_pair(tmp.m_name, std::shared_ptr<DoubleParameter> (new DoubleParameter("mass_"+tmp.m_name,tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) ));
	       gr.insert(std::make_pair(tmp.m_name, std::shared_ptr<DoubleParameter> (new DoubleParameter("width_"+tmp.m_name,tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) ));
	       rr.insert(std::make_pair(tmp.m_name, std::shared_ptr<DoubleParameter> (new DoubleParameter("mag_"+tmp.m_name,tmp.m_strength) ) ));
	       phir.insert(std::make_pair(tmp.m_name, std::shared_ptr<DoubleParameter> (new DoubleParameter("phase_"+tmp.m_name,tmp.m_phase) ) ));

	       unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;
	       //unsigned int last = mr.size()-1;

	       //setup Dynamics
	       //unsigned int last = mr.size()-1;
	       std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
	               *mr[tmp.m_name], *gr[tmp.m_name], tmp.m_mesonRadius, subSys, tmp.m_spin,tmp.m_m,tmp.m_n) );
	       totAmp.addBW(tmpbw, rr[tmp.m_name], phir[tmp.m_name]);

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
		mr.insert(std::make_pair(tmp.m_name, std::shared_ptr<DoubleParameter> (new DoubleParameter("mass_"+tmp.m_name,tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) ));
		gr.insert(std::make_pair(tmp.m_name, std::shared_ptr<DoubleParameter> (new DoubleParameter("width_"+tmp.m_name,tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) ));
		rr.insert(std::make_pair(tmp.m_name, std::shared_ptr<DoubleParameter> (new DoubleParameter("mag_"+tmp.m_name,tmp.m_strength) ) ));
		phir.insert(std::make_pair(tmp.m_name, std::shared_ptr<DoubleParameter> (new DoubleParameter("phase_"+tmp.m_name,tmp.m_phase) ) ));
		DoubleParameter param1("coupling1_"+tmp.m_name,tmp.m_coupling);
		DoubleParameter param2("coupling2_"+tmp.m_name,tmp.m_couplingHidden);

		unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;

		//setup Dynamics
		//unsigned int last = mr.size()-1;
		std::shared_ptr<AmpFlatteRes> tmpbw(new AmpFlatteRes(tmp.m_name.c_str(),
				*mr[tmp.m_name], *gr[tmp.m_name], tmp.m_mesonRadius, param1, param2, \
				PhysConst::instance()->getMass(tmp.m_hiddenParticle1),\
				PhysConst::instance()->getMass(tmp.m_hiddenParticle2),\
				subSys, tmp.m_spin,tmp.m_m,tmp.m_n) );
		totAmp.addBW(tmpbw, rr[tmp.m_name], phir[tmp.m_name]);

		double norm=tmp.m_norm;
		if(norm<0 || _calcNorm) {//recalculate normalization
			norm = normReso(tmpbw);
			reso->m_norm = norm;//updating normalization
		}
		tmpbw->SetNormalization(1/norm);
	}// end loop over resonancesFlatte

	//	ampSetup.save(ampSetup.getFilePath());//save updated information to input file
	nAmps=rr.size();
	if(_calcNorm) integral();
	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity: completed setup!";
}

void AmpSumIntensity::setupTree(allMasses& theMasses, bool isPhspTree){
  DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
  _dpArea = kin->getPhspVolume();
  //  std::cout<<kin->getPhspVolume()<<std::endl;
  _calcNorm=1;
  if(!isPhspTree)
    BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::setupTree() start setting up data Tree";
  else
    BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::setupTree() start setting up phsp Tree";

  //------------Setup Tree---------------------
  std::shared_ptr<FunctionTree> newTree = std::shared_ptr<FunctionTree>(new FunctionTree());

  //------------Setup Tree Pars---------------------
  //std::shared_ptr<DoubleParameter> x = std::shared_ptr<DoubleParameter>(new DoubleParameter("x"));
  std::shared_ptr<MultiDouble> m23 = std::shared_ptr<MultiDouble>( new MultiDouble("m23",theMasses.masses_sq.at( std::make_pair(2,3) )) );
  std::shared_ptr<MultiDouble> m13 = std::shared_ptr<MultiDouble>( new MultiDouble("m13",theMasses.masses_sq.at( std::make_pair(1,3) )) );
  std::shared_ptr<MultiDouble> m12 = std::shared_ptr<MultiDouble>( new MultiDouble("m12",theMasses.masses_sq.at( std::make_pair(1,2) )) );
  treePar = std::shared_ptr<ParameterList>(new ParameterList());
  //treePar->AddParameter(x);
  treePar->AddParameter(m23);
  treePar->AddParameter(m13);
  treePar->AddParameter(m12);

  //----Strategies needed
  std::shared_ptr<MultAll> mmultStrat = std::shared_ptr<MultAll>(new MultAll(ParType::MCOMPLEX));
  std::shared_ptr<AddAll> maddStrat = std::shared_ptr<AddAll>(new AddAll(ParType::MCOMPLEX));
  std::shared_ptr<AbsSquare> msqStrat = std::shared_ptr<AbsSquare>(new AbsSquare(ParType::MDOUBLE));
  std::shared_ptr<LogOf> mlogStrat = std::shared_ptr<LogOf>(new LogOf(ParType::MDOUBLE));
  std::shared_ptr<MultAll> multStrat = std::shared_ptr<MultAll>(new MultAll(ParType::COMPLEX));
  std::shared_ptr<AddAll> addStrat = std::shared_ptr<AddAll>(new AddAll(ParType::DOUBLE));
  std::shared_ptr<AbsSquare> sqStrat = std::shared_ptr<AbsSquare>(new AbsSquare(ParType::DOUBLE));
  std::shared_ptr<LogOf> logStrat = std::shared_ptr<LogOf>(new LogOf(ParType::DOUBLE));
  std::shared_ptr<Complexify> complStrat = std::shared_ptr<Complexify>(new Complexify(ParType::COMPLEX));

  //----Add Top Node and final operations
  newTree->createHead("LH", addStrat); //Sum up all events, collapse multia
  if(!isPhspTree){ //Data: EvtSum of log of Intens needed
    newTree->createNode("Log", mlogStrat, "LH", theMasses.nEvents, false); //log of amp, at each point
    newTree->createNode("Intens", msqStrat, "Log", theMasses.nEvents, false); //I=A^2, at each point
    newTree->createNode("Amplitude", maddStrat, "Intens", theMasses.nEvents, false); //Sum of resonances, at each point
  }else{ //Phsp: PhspSum of Intens needed, log done in LH with other parameters
    newTree->createNode("Intens", msqStrat, "LH", theMasses.nEvents, false); //I=A^2, at each point
    newTree->createNode("Amplitude", maddStrat, "Intens", theMasses.nEvents, false); //Sum of resonances, at each point
  }

  //----Add Resonances
  for(std::vector<Resonance>::iterator reso=ampSetup.getResonances().begin(); reso!=ampSetup.getResonances().end(); reso++){
      Resonance tmp = (*reso);
      if(!tmp.m_enable) continue;
      //setup RooVars
      //namer.push_back(tmp.m_name);
      //mr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("mass_"+tmp.m_name,tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) );
      //gr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("width_"+tmp.m_name,tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) );
     // rr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("rr_"+tmp.m_name,tmp.m_strength) ) );
     // phir.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter("phir_"+tmp.m_name,tmp.m_phase) ) );

      //----Add Nodes
      std::shared_ptr<BreitWignerStrategy> rbwStrat = std::shared_ptr<BreitWignerStrategy>(new BreitWignerStrategy(tmp.m_name,ParType::MCOMPLEX));
      std::shared_ptr<WignerDStrategy> angdStrat = std::shared_ptr<WignerDStrategy>(new WignerDStrategy(tmp.m_name,ParType::MDOUBLE));
      unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;
      unsigned int last = mr.size()-1;
      newTree->createNode("Reso_"+tmp.m_name, mmultStrat, "Amplitude", theMasses.nEvents); //Reso=BW*c*AD
      newTree->createNode("RelBW_"+tmp.m_name, rbwStrat, "Reso_"+tmp.m_name, theMasses.nEvents); //BW
      newTree->createNode("C_"+tmp.m_name, complStrat, "Reso_"+tmp.m_name); //c
      newTree->createLeaf("Intens_"+tmp.m_name, rr[tmp.m_name], "C_"+tmp.m_name); //r
      newTree->createLeaf("Phase_"+tmp.m_name, phir[tmp.m_name], "C_"+tmp.m_name); //phi
      newTree->createNode("AngD_"+tmp.m_name, angdStrat, "Reso_"+tmp.m_name, theMasses.nEvents); //AD
      //BW Par
      newTree->createLeaf("m0_"+tmp.m_name, mr[tmp.m_name], "RelBW_"+tmp.m_name); //m0
      //myTree->createLeaf("x", x, "RelBW_"+tmp.m_name); //x
      switch(subSys){
        case 3:{ //reso in sys of particles 1&2
          //newTree->createLeaf("mym_"+tmp.m_name, m12, "RelBW_"+tmp.m_name); //m
          newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "RelBW_"+tmp.m_name); //ma
          newTree->createLeaf("mb_"+tmp.m_name, kin->m2, "RelBW_"+tmp.m_name); //mb
          break;
        }
        case 4:{ //reso in sys of particles 1&3
          //newTree->createLeaf("mym_"+tmp.m_name, m13, "RelBW_"+tmp.m_name); //m
          newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "RelBW_"+tmp.m_name); //ma
          newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "RelBW_"+tmp.m_name); //mb
          break;
        }
        case 5:{ //reso in sys of particles 2&3
          //newTree->createLeaf("mym_"+tmp.m_name, m23, "RelBW_"+tmp.m_name); //m
          newTree->createLeaf("ma_"+tmp.m_name, kin->m2, "RelBW_"+tmp.m_name); //ma
          newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "RelBW_"+tmp.m_name); //mb
          break;
        }
        default:{
          BOOST_LOG_TRIVIAL(error)<<"AmpSumIntensity::setupTree(): Subsys not found!!";
        }
      }
      newTree->createLeaf("m23", m23, "RelBW_"+tmp.m_name); //ma
      newTree->createLeaf("m13", m13, "RelBW_"+tmp.m_name); //mb
      newTree->createLeaf("m12", m12, "RelBW_"+tmp.m_name); //mc
      newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "RelBW_"+tmp.m_name); //subSysFlag
      newTree->createLeaf("spin_"+tmp.m_name, tmp.m_spin, "RelBW_"+tmp.m_name); //spin
      newTree->createLeaf("d_"+tmp.m_name,  tmp.m_mesonRadius, "RelBW_"+tmp.m_name); //d
      newTree->createLeaf("resWidth_"+tmp.m_name, gr[tmp.m_name], "RelBW_"+tmp.m_name); //resWidth
      newTree->createLeaf("norm_"+tmp.m_name, 1., "RelBW_"+tmp.m_name); //Todo: setup norm, manipulate norm?
      //AD Par
      newTree->createLeaf("m23", m23, "AngD_"+tmp.m_name); //ma
      newTree->createLeaf("m13", m13, "AngD_"+tmp.m_name); //mb
      newTree->createLeaf("m12", m12, "AngD_"+tmp.m_name); //mc
      newTree->createLeaf("M", kin->M, "AngD_"+tmp.m_name); //M
      newTree->createLeaf("m1", kin->m1, "AngD_"+tmp.m_name); //m1
      newTree->createLeaf("m2", kin->m2, "AngD_"+tmp.m_name); //m2
      newTree->createLeaf("m3", kin->m3, "AngD_"+tmp.m_name); //m3
      // unsigned int _subSysFlag = Double_t(paras.GetParameterValue("subSysFlag_"+name));
      newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "AngD_"+tmp.m_name); //subSysFlag
      newTree->createLeaf("spin_"+tmp.m_name,tmp.m_spin, "AngD_"+tmp.m_name); //spin
      newTree->createLeaf("m_"+tmp.m_name, tmp.m_m, "AngD_"+tmp.m_name); //OutSpin 1
      newTree->createLeaf("n_"+tmp.m_name, tmp.m_n, "AngD_"+tmp.m_name); //OutSpin 2

  }// end loop over resonances

  if(isPhspTree) myPhspTree=newTree;
  else myTree=newTree;
}

double AmpSumIntensity::getMaxVal(std::shared_ptr<Generator> gen){
	if(!_calcMaxFcnVal) calcMaxVal(gen);
	return _maxFcnVal;
}

double AmpSumIntensity::getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen){
	calcMaxVal(par,gen);
	return _maxFcnVal;
}

void AmpSumIntensity::calcMaxVal(ParameterList& par, std::shared_ptr<Generator> gen){
	setParameterList(par);
	return calcMaxVal(gen);
}

void AmpSumIntensity::calcMaxVal(std::shared_ptr<Generator> gen){
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	//	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	//	boost::minstd_rand rndGen2(123);
	//	boost::uniform_real<> uni_dist13(xLimit_low[0],xLimit_high[0]);
	//	boost::uniform_real<> uni_dist23(xLimit_low[1],xLimit_high[1]);
	//	boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uni13(rndGen2, uni_dist13);
	//	boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uni23(rndGen2, uni_dist23);
	double maxM23=-999; double maxM13=-999; double maxVal=0;
	for(unsigned int i=0; i<100000; i++){
		double m23sq=gen->getUniform()*(kin->m23_sq_max-kin->m23_sq_min)+kin->m23_sq_min;
		double m13sq=gen->getUniform()*(kin->m13_sq_max-kin->m13_sq_min)+kin->m13_sq_min;
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
			<<amp->GetName()<<": "<<1.0/norm;
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

std::complex<double> AmpSumIntensity::getFirstBW(dataPoint& point, ParameterList& par){
    setParameterList(par);
    return totAmp.getFirstBW(point);
}

std::complex<double> AmpSumIntensity::getFirstReso(dataPoint& point, ParameterList& par){
    setParameterList(par);
    return totAmp.getFirstReso(point);
}

std::complex<double> AmpSumIntensity::getFirstAmp(dataPoint& point, ParameterList& par){
    setParameterList(par);
    return totAmp.getFirstAmp(point);
}

const ParameterList& AmpSumIntensity::intensity(std::vector<double> point, ParameterList& par){
	setParameterList(par);
	//	dataPoint dataP; dataP.setVal("m23sq",point[0]); dataP.setVal("m13sq",point[1]);
	dataPoint dataP; dataP.setVal(0,point[0]); dataP.setVal(1,point[1]);
	return intensity(dataP);
}
const ParameterList& AmpSumIntensity::intensity(dataPoint& point, ParameterList& par){
	setParameterList(par);
	return intensity(point);
}
const ParameterList& AmpSumIntensity::intensity(dataPoint& point){
	double AMPpdf=0;
	if(Kinematics::instance()->isWithinPhsp(point)) AMPpdf = totAmp.evaluate(point);
//	else BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity: data point "
//			<<" (m13sq="<<point.getVal(4)
//			<<" m23sq="<<point.getVal(5)<<")"
//			<<" not within phsp! Skip!";//shall we give an error if datapoint is outside phsp?

	if(AMPpdf!=AMPpdf){
		BOOST_LOG_TRIVIAL(error)<<"Error AmpSumIntensity: Intensity is not a number!!";
		AMPpdf = 0;
	}
	double eff=eff_->evaluate(point);
	result.SetParameterValue(0,AMPpdf*eff);
	return result;
}
std::shared_ptr<FunctionTree> AmpSumIntensity::functionTree(allMasses& theMasses) {
	/*if(outPar.GetNParameter()>0) return std::shared_ptr<FunctionTree>();
	fillStartParVec(outPar);
	//outPar.AddParameter(treePar->GetDoubleParameter("x"));
	outPar.AddParameter(treePar->GetDoubleParameter("m23"));
	outPar.AddParameter(treePar->GetDoubleParameter("m13"));
	outPar.AddParameter(treePar->GetDoubleParameter("m12"));*/

	if(myTree) return myTree;
	else setupTree(theMasses, false);

	return myTree;
}
std::shared_ptr<FunctionTree> AmpSumIntensity::phspTree(allMasses& theMasses) {
    /*if(outPar.GetNParameter()>0) return std::shared_ptr<FunctionTree>();
    fillStartParVec(outPar);
    //outPar.AddParameter(treePar->GetDoubleParameter("x"));
    outPar.AddParameter(treePar->GetDoubleParameter("m23"));
    outPar.AddParameter(treePar->GetDoubleParameter("m13"));
    outPar.AddParameter(treePar->GetDoubleParameter("m12"));*/

    if(myPhspTree) return myPhspTree;
    else setupTree(theMasses, true);

    return myPhspTree;
}
void AmpSumIntensity::setParameterList(ParameterList& par){
	//parameters varied by Minimization algorithm
	for(unsigned int i=0; i<nAmps; i++){
//		*rr[i] = DoubleParameter(par.GetDoubleParameter(2*i));
//		*phir[i] = DoubleParameter(par.GetDoubleParameter(2*i+1));
	    if(!rr[namer[i]]->IsFixed()){
	      rr[namer[i]]->SetValue(par.GetDoubleParameter(2*i)->GetValue());
		  rr[namer[i]]->SetError(par.GetDoubleParameter(2*i)->GetError());
	    }
	    if(!phir[namer[i]]->IsFixed()){
		  phir[namer[i]]->SetValue(par.GetDoubleParameter(2*i+1)->GetValue());
		  phir[namer[i]]->SetError(par.GetDoubleParameter(2*i+1)->GetError());
	    }
	}
	return;
}
const bool AmpSumIntensity::fillStartParVec(ParameterList& outPar){
	if(outPar.GetNParameter())
		return false; //already filled, TODO: exception?
	for(unsigned int i=0; i<namer.size();i++){
		//add strength and phases of the used amplitudes
		outPar.AddParameter(rr[namer[i]]);
		outPar.AddParameter(phir[namer[i]]);
	}
	return true;
}

void AmpSumIntensity::printAmps(){
	std::stringstream outStr;
	outStr<<"AmpSumIntensity: Printing amplitudes with current(!) set of parameters:\n";
	for(unsigned int i=0;i<nAmps;i++){
		outStr<<std::setw(12)<<totAmp.getAmpName(i)<<": ";
		outStr<<"Mag = "<<std::setw(6)<<rr[namer[i]]->GetValue()<<" +- "<<std::setw(5)<<rr[namer[i]]->GetError()->GetError();
		outStr<<"   Phase = "<<std::setw(6)<<phir[namer[i]]->GetValue()<<" +- "<<std::setw(5)<<phir[namer[i]]->GetError()->GetError();
		if(!(i==nAmps-1)) outStr << "\n";
	}
	BOOST_LOG_TRIVIAL(info)<<outStr.str();
	return;
}
void AmpSumIntensity::printFractions(){
	std::stringstream outStr;
	outStr<<"Fit fractions for all amplitudes: \n";
	double sumFrac=0;
	for(unsigned int i=0;i<totAmp.getNAmps();i++){
		double frac = getFraction(i);
		sumFrac+=frac;
		outStr<<std::setw(10)<<totAmp.getAmpName(i)<<":    "<<frac<<"\n";
		//		if(!(i==totAmp.getNAmps()-1)) outStr << "\n";
	}
	outStr<<std::setw(10)<<" "<<"    ==========\n";
	outStr<<std::setw(10)<<" "<<"     "<<sumFrac;
	BOOST_LOG_TRIVIAL(info)<<outStr.str();
	return;
}

double AmpSumIntensity::getIntValue(std::string var1, double min1, double max1, std::string var2, double min2, double max2){
	/*
	 * Integrates in \var1 from \min1 to \max1 and in \var2 from \min2 to \max2.
	 * Is intended to be used for calculation of bin content.
	 */
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double _min1 = min1;
	double _min2 = min2;
	double _max1 = max1;
	double _max2 = max2;
	//if(_min1==0) _min1 = kin->GetMin(var1);
	//if(_max1==0) _min1 = kin->GetMax(var1);
	if(_min2==0) _min2 = kin->getMin(var2);
	if(_max2==0) _max2 = kin->getMax(var2);
	unsigned int dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	double xLimit_low[2];
	double xLimit_high[2];

	if(var1 == "m13sq" && var2 == "m23sq"){
		xLimit_low[0] = _min1;
		xLimit_low[1] = _min2;
		xLimit_high[0] = _max1;
		xLimit_high[1] = _max2;
	}else if(var1 == "m23sq" && var2 == "m13sq"){
		xLimit_low[0] = _min2;
		xLimit_low[1] = _min1;
		xLimit_high[0] = _max2;
		xLimit_high[1] = _max1;
	} else {
		BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::getIntValue() wrong variables specified!";
		return -999;
	}

	size_t calls = 5000;
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

	return res;
}
