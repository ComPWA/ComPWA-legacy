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

      //setup Dynamics and Angular Distribution
      unsigned int last = mr.size()-1;
      if(tmp.m_daugtherA==2 && tmp.m_daugtherB==3){
        std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
            tmp.m_name.c_str(), ma, *mr[last], *gr[last], *qr[last], 1, tmp.m_spin) );
        tmpbw->setDecayMasses(m2, m3);
        rbw.push_back(tmpbw);
        angd.push_back( std::shared_ptr<AmpWigner> (new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
            mc, ma, mb, 1, *aj[last], *am[last], *an[last]) ) );
        totAmp.addBW(rbw.at(last), *rr.at(last), *phir.at(last), angd.at(last));
      }else if(tmp.m_daugtherA==1 && tmp.m_daugtherB==3){
        std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
            tmp.m_name.c_str(), mb, *mr[last], *gr[last], *qr[last], 2, tmp.m_spin) );
        tmpbw->setDecayMasses(m1, m3);
        rbw.push_back(tmpbw);
        angd.push_back( std::shared_ptr<AmpWigner> (new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
            mc, ma, mb, 2, *aj[last], *am[last], *an[last]) ) );
        totAmp.addBW(rbw.at(last), *rr.at(last), *phir.at(last), angd.at(last));
      }else if(tmp.m_daugtherA==1 && tmp.m_daugtherB==2){
        std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
            tmp.m_name.c_str(), mc, *mr[last], *gr[last], *qr[last], 3, tmp.m_spin) );
        tmpbw->setDecayMasses(m1, m2);
        rbw.push_back(tmpbw);
        angd.push_back( std::shared_ptr<AmpWigner> (new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
            mc, ma, mb, 3, *aj[last], *am[last], *an[last]) ) );
        totAmp.addBW(rbw.at(last), *rr.at(last), *phir.at(last), angd.at(last));
      }else{ //ignore resonance
          //std::cout << "Problem" << std::cout;
          mr.pop_back();
          qr.pop_back();
          gr.pop_back();
          rr.pop_back();
          phir.pop_back();
          aj.pop_back();
          am.pop_back();
          an.pop_back();
      }

      //rbw.at(last).printName(std::cout);
      //angd.at(last).printName(std::cout);
    }

    std::cout << "completed setup" << std::endl;
  }

  virtual const double integral(ParameterList& par){
    double integral = 0;
    /*unsigned int nSteps = 1000000;   TODO
    double step = (max_-min_)/(double)nSteps;

    integral += step*BreitWigner(min_, par.at(0)->GetValue(), par.at(1)->GetValue())/2.;
    for(unsigned int k=1; k<nSteps; k++){
      integral += step*BreitWigner((min_+k*step), par.at(0)->GetValue(), par.at(1)->GetValue());
    }
    integral += step*BreitWigner(max_, par.at(0)->GetValue(), par.at(1)->GetValue())/2.;*/

    return integral;
  }

  virtual const double intensity(std::vector<double> x, ParameterList& par){
    ma.setVal(sqrt(x[0])); mb.setVal(sqrt(x[1])); mc.setVal(sqrt(x[2]));
    double AMPpdf = totAmp.evaluate();
    if(AMPpdf!=AMPpdf){
      //cout << "Error: Intensity is not a number!!" << endl; TODO: exception
      return 0;
    }
    return AMPpdf;
  }

  virtual const bool fillStartParVec(ParameterList& outPar){
    if(outPar.GetNParameter())
      return false; //already filled, TODO: exception?

    outPar.AddParameter(DoubleParameter(1.5, 0.5, 2.5, 0.1));
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
  std::vector<std::shared_ptr<RooRealVar> > phir;
  std::vector<std::shared_ptr<RooRealVar> > aj;
  std::vector<std::shared_ptr<RooRealVar> > am;
  std::vector<std::shared_ptr<RooRealVar> > an;

  std::vector<std::shared_ptr<AmpAbsDynamicalFunction> > rbw;
  std::vector<std::shared_ptr<AmpWigner> > angd;

  /*AmpRelBreitWignerRes * rbwm23_0;
  AmpRelBreitWignerRes * rbwm23_1;
  AmpRelBreitWignerRes * rbwm23_2;
  AmpRelBreitWignerRes * inter13_0;
  AmpRelBreitWignerRes * inter12_0;*/

 /* AmpWigner * angd23_0;
  AmpWigner * angd23_1;
  AmpWigner * angd13_i;
  AmpWigner * angd12_i;*/

 private:


};

#endif
