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

#include "Amplitude.hpp"
#include "PWAParameter.hpp"
#include "PWAGenericPar.hpp"

#include "AmpRelBreitWignerRes.hpp"
#include "AmpGausRes.hpp"
#include "AmpFlatteRes.hpp"
#include "AmpWigner.hpp"
#include "AmpSumOfAmplitudes.hpp"

class AmpSumIntensity : public Amplitude {

public:
  /// Default Constructor (0x0)
  AmpSumIntensity(const double inM, const double inBr, const double in1, const double in2, const double in3)
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
    totAmp23("rbwAmp23", "totAmp"){

    //AmpSumOfAmplitudes totAmp23("rbwAmp23", "totAmp");

    //Resonance Parameter
    mr0 = RooRealVar("m_{f0(1710)}", "resMass", 1.5, 1.0, 2.0);
    qr0 = RooRealVar("q0", "breakupMom", 1.0);
    gr0 = RooRealVar("#Gamma_{f(1710)}", "resWidth", 0.05, 0., 2.) ;
    rr0 = RooRealVar("r_{0}", "amplitude", 1.);
    phir0 = RooRealVar("#phi_{0}", "phase", 2.*PI/3.);

    mr1 = RooRealVar("m_{fJ(2220)}", "resMass", 2.1, 1.5, 2.9);
    qr1 = RooRealVar("q1", "breakupMom", 1.0);
    gr1 = RooRealVar("#Gamma_{f(2220)}", "resWidth", 0.1, 0., 2.) ;
    rr1 = RooRealVar("r_{1}", "amplitude", 1);
    phir1 = RooRealVar("#phi_{1}", "phase", PI);

    mr2 = RooRealVar("m_{fJ(1900)}", "resMass", 2.0, 1.5, 2.5);
    qr2 = RooRealVar("q2", "breakupMom", 1.0);
    gr2 = RooRealVar("#Gamma_{f(1900)}", "resWidth", 0.1, 0., 2.) ;
    rr2 = RooRealVar("r_{2}", "amplitude", 1);
    phir2 = RooRealVar("#phi_{2}", "phase", 2*PI);

    mi0 = RooRealVar("m_{Hyp(1500)}", "resMass", 1., 0.8, 1.2);
    qi0 = RooRealVar("qi", "breakupMom", 1.0);
    gi0 = RooRealVar("#Gamma_{Hyp(1500)}", "resWidth", 0.05, 0., 8.) ;
    ri0 = RooRealVar("r_{i}", "amplitude", 0.5);
    phii0 = RooRealVar("#phi_{i}", "phase", 0);

    //Resonances
    rbwm23_0 = new AmpRelBreitWignerRes("f1710", "f1710", ma, mr0, gr0, qr0, 1, 0);
    //AmpFlatteRes * rbwm23_0 = new AmpFlatteRes("f1710", "f1710", ma, mr0, gr0, qr0, 1, 0);
    rbwm23_1 = new AmpRelBreitWignerRes("f2220", "f2220", ma, mr1, gr1, qr1, 1, 0);
    rbwm23_2 = new AmpRelBreitWignerRes("f1900", "f1900", ma, mr2, gr2, qr2, 1, 0);
    //AmpGausRes * rbwm23_1 = new AmpGausRes("f2220", "f2220", ma, mr1, gr1, 0);

    inter13_0 = new AmpRelBreitWignerRes("H1500_a", "H1500_a", mb, mi0, gi0, qi0, 2, 0);

    inter12_0 = new AmpRelBreitWignerRes("H1500_b", "H1500_b", mc, mi0, gi0, qi0, 3, 0);

    //const double mPi = 0.13957;
    rbwm23_0->setDecayMasses (m2, m3);
  //  rbwm23_0->setBarrierMass (0.8, 0.8);
    rbwm23_1->setDecayMasses (m2, m3);
    rbwm23_2->setDecayMasses (m2, m3);

    inter13_0->setDecayMasses (m1, m3);

    inter12_0->setDecayMasses (m1, m2);

    //Angular Distribution Parameter
    aj0 = RooRealVar("j_{f(1710)}", "j", 2.);
    am0 = RooRealVar("m_{f(1710)}", "m", 0.);
    an0 = RooRealVar("n_{f(1710)}", "n", 0.) ;

    aj1 = RooRealVar("j_{f(2220)}", "j", 1.);
    am1 = RooRealVar("m_{f(2220)}", "m", 0.);
    an1 = RooRealVar("n_{f(2220)}", "n", 0.);

    aji = RooRealVar("j_{H(1500)}", "j", 0.);
    ami = RooRealVar("m_{H(1500)}", "m", 0.);
    ani = RooRealVar("n_{H(1500)}", "n", 0.);

    //Angular Distributions
    angd23_0 = new AmpWigner("a_f1710", "a_f1710", mc, ma, mb, 1, aj0, am0, an0);
    angd23_1 = new AmpWigner("a_f2220", "a_f2220", mc, ma, mb, 1, aj1, am1, an1);

    angd13_i = new AmpWigner("a_H1500_a", "a_H1500_a", mc, ma, mb, 2, aji, ami, ani);

    angd12_i = new AmpWigner("a_H1500_b", "a_H1500_b", mc, ma, mb, 3, aji, ami, ani);

    totAmp23.addBW (rbwm23_0, rr0, phir0, angd23_0);
    totAmp23.addBW (rbwm23_1, rr1, phir1, angd23_0);
    totAmp23.addBW (rbwm23_2, rr2, phir2, angd23_1);
    totAmp23.addBW (inter13_0, ri0, phii0, angd13_i);
    totAmp23.addBW (inter12_0, ri0, phii0, angd12_i);
  }

  virtual const double integral(std::vector<std::shared_ptr<PWAParameter> >& par){
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

  virtual const double intensity(std::vector<double> x, std::vector<std::shared_ptr<PWAParameter> >& par){
    ma.setVal(sqrt(x[0])); mb.setVal(sqrt(x[1])); mc.setVal(sqrt(x[2]));
    double AMPpdf = totAmp23.evaluate();
    if(AMPpdf!=AMPpdf){
      //cout << "Error: Intensity is not a number!!" << endl;
      return 0;
    }
    return AMPpdf;
  }
  virtual const bool fillStartParVec(std::vector<std::shared_ptr<PWAParameter> >& outPar){
    if(outPar.size())
      return false; //already filled
    outPar.push_back(shared_ptr<PWAParameter>(new PWAGenericPar<double>(1.5, 0.5, 2.5, 0.1)));
    return true;
  }

  /** Destructor */
  virtual ~AmpSumIntensity(){
    delete rbwm23_0;
    delete rbwm23_1;
    delete rbwm23_2;
    delete inter13_0;
    delete inter12_0;
    delete angd23_0;
    delete angd23_1;
    delete angd13_i;
    delete angd12_i;
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

  AmpSumOfAmplitudes totAmp23;

  //Resonance Variables
  RooRealVar mr0;
  RooRealVar qr0;
  RooRealVar gr0;
  RooRealVar rr0;
  RooRealVar phir0;

  RooRealVar mr1;
  RooRealVar qr1;
  RooRealVar gr1;
  RooRealVar rr1;
  RooRealVar phir1;

  RooRealVar mr2;
  RooRealVar qr2;
  RooRealVar gr2;
  RooRealVar rr2;
  RooRealVar phir2;

  RooRealVar mi0;
  RooRealVar qi0;
  RooRealVar gi0;
  RooRealVar ri0;
  RooRealVar phii0;

  AmpRelBreitWignerRes * rbwm23_0;
  AmpRelBreitWignerRes * rbwm23_1;
  AmpRelBreitWignerRes * rbwm23_2;
  AmpRelBreitWignerRes * inter13_0;
  AmpRelBreitWignerRes * inter12_0;

  //Angular Distribution Variables
  RooRealVar aj0;
  RooRealVar am0;
  RooRealVar an0;

  RooRealVar aj1;
  RooRealVar am1;
  RooRealVar an1;

  RooRealVar aji;
  RooRealVar ami;
  RooRealVar ani;

  AmpWigner * angd23_0;
  AmpWigner * angd23_1;
  AmpWigner * angd13_i;
  AmpWigner * angd12_i;

 private:


};

#endif
