//****************************************************************************
// Wrapper to provide intensity of amplitude sum
//****************************************************************************

#ifndef _AMPSUMINTENSITY_HPP
#define _AMPSUMINTENSITY_HPP

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

#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"

class AmpSumIntensity : public Amplitude {

public:
  /// Default Constructor (0x0)
  AmpSumIntensity(const double inM, const double inBr, const double in1
      ,const double in2, const double in3, AmplitudeSetup ini)
  : M(inM), Br(inBr), m1(in1), m2(in2), m3(in3),PI(3.14159),
    m23_sq_min((m2+m3)*(m2+m3)),
    m23_sq_max((M-m1)*(M-m1)),
    m13_sq_min((m1+m3)*(m1+m3)),
    m13_sq_max((M-m2)*(M-m2)),
    m12_sq_min((m1+m2)*(m1+m2)),
    m12_sq_max((M-m3)*(M-m3)),
    m23_min((m2+m3)), m23_max((M-m1)),
    m13_min((m1+m3)), m13_max((M-m2)),
    m12_min((m1+m2)), m12_max((M-m3)),
    ma("ma", "mass", m23_min, m23_max),
    mb("mb", "mass", m13_min, m13_max),
    mc("mc", "mass", m12_min, m12_max),
    totAmp("relBWsumAmplitude", "totAmp"){

    //AmpSumOfAmplitudes totAmp23("rbwAmp23", "totAmp");


    //unsigned int numReso = ini.getResonances().size();
    for(std::vector<Resonance>::iterator reso=ini.getResonances().begin(); reso!=ini.getResonances().end(); reso++){
      Resonance tmp = (*reso);
      //setup RooVars
      mr.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("m_{"+tmp.m_name+"}").c_str(), "resMass", tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) );
      qr.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("q_{"+tmp.m_name+"}").c_str(), "breakupMom", tmp.m_breakup_mom) ) );
      gr.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("g_{"+tmp.m_name+"}").c_str(), "resWidth", tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) );
      rr.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("r_{"+tmp.m_name+"}").c_str(), "amplitude", tmp.m_strength) ) );
      phir.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("phi_{"+tmp.m_name+"}").c_str(), "phase", tmp.m_phase) ) );

      aj.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("j_{"+tmp.m_name+"}").c_str(), "j", tmp.m_spin) ) );
      am.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("m_{"+tmp.m_name+"}").c_str(), "m", tmp.m_m) ) );
      an.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("n_{"+tmp.m_name+"}").c_str(), "n", tmp.m_n) ) );


//      unsigned int last = mr.size()-1;
//      if(tmp.m_daugtherA==2 && tmp.m_daugtherB==3){
//        std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
//            tmp.m_name.c_str(), ma, *mr[last], *gr[last], *qr[last], 1, tmp.m_spin) );
//        tmpbw->setDecayMasses(m2, m3);
//        rbw.push_back(tmpbw);
//        angd.push_back( std::shared_ptr<AmpWigner> (new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
//            mc, ma, mb, 1, *aj[last], *am[last], *an[last]) ) );
//        totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last), angd.at(last));
//      }else if(tmp.m_daugtherA==1 && tmp.m_daugtherB==3){
//        std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
//            tmp.m_name.c_str(), mb, *mr[last], *gr[last], *qr[last], 2, tmp.m_spin) );
//        tmpbw->setDecayMasses(m1, m3);
//        rbw.push_back(tmpbw);
//        angd.push_back( std::shared_ptr<AmpWigner> (new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
//            mc, ma, mb, 2, *aj[last], *am[last], *an[last]) ) );
//        totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last), angd.at(last));
//      }else if(tmp.m_daugtherA==1 && tmp.m_daugtherB==2){
//        std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
//            tmp.m_name.c_str(), mc, *mr[last], *gr[last], *qr[last], 3, tmp.m_spin) );
//        tmpbw->setDecayMasses(m1, m2);
//        rbw.push_back(tmpbw);
//        angd.push_back( std::shared_ptr<AmpWigner> (new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
//            mc, ma, mb, 3, *aj[last], *am[last], *an[last]) ) );
//        totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last), angd.at(last));
//      }else{ //ignore resonance
//          //std::cout << "Problem" << std::cout;
//          mr.pop_back();
//          qr.pop_back();
//          gr.pop_back();
//          rr.pop_back();
//          phir.pop_back();
//          aj.pop_back();
//          am.pop_back();
//          an.pop_back();
//      }
      par1.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("par1_{"+tmp.m_name+"}").c_str(), "n", tmp.m_par1) ) );
      par2.push_back( std::shared_ptr<RooRealVar> (new RooRealVar(("par2_{"+tmp.m_name+"}").c_str(), "n", tmp.m_par2) ) );
      const double* mass_first;
      const double* mass_third;
      const double* mass_second;
      RooAbsReal* mass_eval_12;
      RooAbsReal* mass_eval_23;
      RooAbsReal* mass_eval_13;
      int count;
      if(tmp.m_daugtherA==2 && tmp.m_daugtherB==3){
          mass_first = &m2;
          mass_second = &m3;
          mass_third = &m1;
          mass_eval_12 = &ma;
          mass_eval_23 = &mb;
          mass_eval_13 = &mc;
          count=1;
      }else if(tmp.m_daugtherA==1 && tmp.m_daugtherB==3){
          mass_first = &m1;
          mass_second = &m3;
          mass_third = &m2;
          mass_eval_12 = &mb;
          mass_eval_23 = &mc;
          mass_eval_13 = &ma;
          count=2;
      }else if(tmp.m_daugtherA==1 && tmp.m_daugtherB==2){
          mass_first = &m1;
          mass_second = &m2;
          mass_third = &m3;
          mass_eval_12 = &mc;
          mass_eval_23 = &ma;
          mass_eval_13 = &mb;
          count=3;
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
          std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
                  tmp.m_name.c_str(), *mass_eval_12, *mr[last], *gr[last], *qr[last], 1, tmp.m_spin) );
    	  tmpbw->setDecayMasses(*mass_first,*mass_second);
          rbw.push_back(tmpbw);
//         angd.push_back( std::shared_ptr<AmpWigner> (new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
//			  *mass_eval_23, *mass_eval_12, *mass_eval_13, count, *aj[last], *am[last], *an[last]) ) );
          std::shared_ptr<AmpWigner> tmpWigner(new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
			  mc, ma, mb, count, *aj[last], *am[last], *an[last]) ) ;
          tmpWigner->setDecayMasses(M, *mass_first, *mass_second , *mass_third);
          angd.push_back( tmpWigner );
          totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last), angd.at(last));
//          totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last));
      }
      else if(tmp.m_type=="flatte"){
          std::shared_ptr<AmpFlatteRes> tmpbw(new AmpFlatteRes(tmp.m_name.c_str(),
                  tmp.m_name.c_str(), *mass_eval_12, *mr[last], *gr[last], *qr[last], *par1[last], *par2[last], count, tmp.m_spin) );
          tmpbw->setDecayMasses(*mass_first,*mass_second);
          tmpbw->setBarrierMass(0.547853,0.1396);//a_0->eta pi hidden channel
          rbw.push_back(tmpbw);
          std::shared_ptr<AmpWigner> tmpWigner(new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
			  mc, ma, mb, count, *aj[last], *am[last], *an[last]) ) ;
          tmpWigner->setDecayMasses(M, *mass_first, *mass_second , *mass_third);
          angd.push_back( tmpWigner );
          totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last), angd.at(last));
//    	  totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last));
          cout<<"12123"<<endl;
      }
      else continue;
    }

    std::cout << "completed setup" << std::endl;
  }

  virtual const double integral(ParameterList& par){
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

  virtual const double intensity(std::vector<double>& x, ParameterList& par){
    //TODO: check x exception
    if(x.size()!=3) return 0;

    ma.setVal(x[0]); mb.setVal(x[1]); mc.setVal(x[2]);

    if( par.GetNDouble()>rr.size() ){
        std::cout << "Error: Parameterlist doesn't match model!!" << std::endl; //TODO: exception
      return 0;
    }

	//parameters varied by Minimization algorithm
    for(unsigned int i=0; i<rr.size(); i++){
      rr[i]->setVal(par.GetDoubleParameter(i).GetValue());//free
      //phir[i]->setVal(par.GetDoubleParameter(2*i+1).GetValue());//fixed
    }
    //rr[rr.size()-1]->setVal(par.GetDoubleParameter(2*(rr.size()-1)).GetValue());

    double AMPpdf = totAmp.evaluate();

    if(AMPpdf!=AMPpdf){
        //std::cout << "Error: Intensity is not a number!!" << std::endl; //TODO: exception
      return 0;
    }

    return AMPpdf;
  }

  virtual const bool fillStartParVec(ParameterList& outPar){
    if(outPar.GetNParameter())
      return false; //already filled, TODO: exception?

    for(unsigned int i=0; i<rr.size();i++){
      //add strength and phases of the used amplitudes

      //outPar.AddParameter(DoubleParameter(rr[i]->getVal()));
      if(rr[i]->hasError()) //TODO: check bounds
        outPar.AddParameter(DoubleParameter(rr[i]->getVal(), rr[i]->getMin(), rr[i]->getMax(), rr[i]->getError()));
      else
        outPar.AddParameter(DoubleParameter(rr[i]->getVal(), 0.1));
      //outPar.AddParameter(DoubleParameter(phir[i]->getVal(), phir[i]->getMin(), phir[i]->getMax(), phir[i]->getError()));
    }
    //outPar.AddParameter(DoubleParameter(rr[rr.size()-1]->getVal(), rr[rr.size()-1]->getMin(), rr[rr.size()-1]->getMax(), rr[rr.size()-1]->getError()));
    //outPar.AddParameter(DoubleParameter(1.5, 0.5, 2.5, 0.1));
    return true;
  }

  /** Destructor */
  virtual ~AmpSumIntensity(){

  }

 protected:
//constants
  const Double_t M; // GeV/c² (J/psi+)
  const Double_t Br; // GeV/c² (width)
  const Double_t m1; // GeV/c² (gamma)
  const Double_t m2; // GeV/c² (pi)
  const Double_t m3; // GeV/c² (pi)
  //const Double_t c = 299792458.; // m/s
  const Double_t PI; // m/s

  const Double_t m23_sq_min;
  const Double_t m23_sq_max;
  const Double_t m13_sq_min;
  const Double_t m13_sq_max;
  const Double_t m12_sq_min;
  const Double_t m12_sq_max;

  const Double_t m23_min;
  const Double_t m23_max;
  const Double_t m13_min;
  const Double_t m13_max;
  const Double_t m12_min;
  const Double_t m12_max;

  RooRealVar ma;
  RooRealVar mb;
  RooRealVar mc;

  AmpSumOfAmplitudes totAmp;

  //Resonance Variables
  std::vector<std::shared_ptr<RooRealVar> > mr;
  std::vector<std::shared_ptr<RooRealVar> > qr;
  std::vector<std::shared_ptr<RooRealVar> > gr;
  std::vector<std::shared_ptr<RooRealVar> > rr;
  std::vector<std::shared_ptr<RooRealVar> > par1;
  std::vector<std::shared_ptr<RooRealVar> > par2;
  std::vector<std::shared_ptr<RooRealVar> > phir;
  std::vector<std::shared_ptr<RooRealVar> > aj;
  std::vector<std::shared_ptr<RooRealVar> > am;
  std::vector<std::shared_ptr<RooRealVar> > an;

  std::vector<std::shared_ptr<AmpAbsDynamicalFunction> > rbw;
  std::vector<std::shared_ptr<AmpWigner> > angd;

  double lambda(double x, double y, double z){
    return x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z;
  }

  double m13_sq_max_constr(double &m23_sq){
    return m1*m1+m3*m3+0.5/m23_sq*((M*M-m23_sq-m1*m1)*(m23_sq-m2*m2+m3*m3)+sqrt(lambda(m23_sq,M*M,m1*m1))*sqrt(lambda(m23_sq,m2*m2,m3*m3)));
  }

  double m13_sq_min_constr(double &m23_sq){
    return m1*m1+m3*m3+0.5/m23_sq*((M*M-m23_sq-m1*m1)*(m23_sq-m2*m2+m3*m3)-sqrt(lambda(m23_sq,M*M,m1*m1))*sqrt(lambda(m23_sq,m2*m2,m3*m3)));
  }

  double m23_sq_max_constr(double &m13_sq){
    return m2*m2+m3*m3+0.5/m13_sq*((M*M-m13_sq-m2*m2)*(m13_sq-m1*m1+m3*m3)+sqrt(lambda(m13_sq,M*M,m2*m2))*sqrt(lambda(m13_sq,m1*m1,m3*m3)));
  }

  double m23_sq_min_constr(double &m13_sq){
    return m2*m2+m3*m3+0.5/m13_sq*((M*M-m13_sq-m2*m2)*(m13_sq-m1*m1+m3*m3)-sqrt(lambda(m13_sq,M*M,m2*m2))*sqrt(lambda(m13_sq,m1*m1,m3*m3)));
  }

  double mb_sq_max_constr(double &ma_sq){
    return m13_sq_max_constr(ma_sq);
  }

  double mb_sq_min_constr(double &ma_sq){
    return m13_sq_min_constr(ma_sq);
  }

  double ma_sq_max_constr(double &mb_sq){
    return m23_sq_max_constr(mb_sq);
  }

  double ma_sq_min_constr(double &mb_sq){
    return m23_sq_min_constr(mb_sq);
  }


 private:


};

#endif
