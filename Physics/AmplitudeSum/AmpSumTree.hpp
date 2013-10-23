//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//****************************************************************************
// Wrapper to provide intensity of amplitude sum
//****************************************************************************

#ifndef _AMPSUMTREE_HPP
#define _AMPSUMTREE_HPP

#include <vector>
#include <memory>

#include "Physics/Amplitude.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Functions.hpp"

#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"

class BreitWignerStrategy : public Strategy {
public:
  BreitWignerStrategy(const std::string resonanceName):name(resonanceName){
    //name = +resonanceName;
  }

  virtual const std::string to_str() const {
    return ("relativistic BreitWigner of "+name);
  }

  virtual std::shared_ptr<AbsParameter> execute(const ParameterList& paras) {

    double Gamma0, GammaV;
    double m0 = Double_t(paras.GetParameterValue("m0_"+name));
    double m  = Double_t(paras.GetParameterValue("x_"+name));
    double ma = Double_t(paras.GetParameterValue("ma_"+name));
    double mb = Double_t(paras.GetParameterValue("mb_"+name));
    unsigned int spin = Double_t(paras.GetParameterValue("spin_"+name));
    double d = Double_t(paras.GetParameterValue("d_"+name));

    Gamma0 = double(paras.GetParameterValue("resWidth_"+name));
    GammaV = Gamma0 * (m0 / m) * pow(q(ma,mb,m) / q0(ma,mb,m0), 2.*spin + 1.)  * BLprime2(ma,mb,m0,m,d,spin);

    RooComplex denom = RooComplex(m0*m0 - m*m, -m0 * GammaV);
    RooComplex res = RooComplex(m0 * Gamma0) / denom;

    std::complex<double> result (res.re(),res.im());
    std::shared_ptr<ComplexParameter> bw(new ComplexParameter("relBW of "+name, result));
    return bw;
  }

protected:
  std::string name;

  double q0(const double& ma, const double& mb, const double& m0) const {
    double mapb = ma + mb;
    double mamb = ma - mb;

    return sqrt ( (m0*m0 - mapb*mapb) * (m0*m0 - mamb*mamb) ) / (2. * m0 );
  }

  double q(const double& ma, const double& mb, const double& x) const {
    double mapb = ma + mb;
    double mamb = ma - mb;

    return sqrt ( (x*x - mapb*mapb) * (x*x - mamb*mamb) ) / (2. * x );
  }


  // compute part of the Blatt-Weisskopf barrier factor
  //   BLprime = sqrt (F(q0)/F(q))
  double F(const double& p, const double& d, unsigned int& spin) const {
    double retVal = 1;

    if (spin == 0)
      retVal = 1;
    else if (spin == 1)
      retVal = 1 + p*p * d*d;
    else if (spin == 2) {
      double z = p*p * d*d;
      retVal = (z-3.)*(z-3.) + 9*z;
    }
    return retVal;
  }


  // compute square of Blatt-Weisskopf barrier factor
  double BLprime2(const double& ma, const double& mb, const double& m0, const double& x, const double& d, unsigned int& spin) const {
    //  cout << q0() << " " << q() << "\t" << F(q0()) << " " << F(q()) << endl;
    return F(q0(ma, mb, m0),d,spin) / F(q(ma, mb, x),d,spin);
  }

};

class WignerDStrategy : public Strategy {
public:
  WignerDStrategy(const std::string resonanceName):name(resonanceName){
    //name = +resonanceName;
  }

  virtual const std::string to_str() const {
    return ("WignerD of "+name);
  }

  virtual std::shared_ptr<AbsParameter> execute(const ParameterList& paras) {

    double Gamma0, GammaV;
    double _m23 = Double_t(paras.GetParameterValue("m23"));
    double _m13 = Double_t(paras.GetParameterValue("m13"));
    double _m12 = Double_t(paras.GetParameterValue("m12"));
    double _M  = Double_t(paras.GetParameterValue("M"));
    double _m1 = Double_t(paras.GetParameterValue("m1"));
    double _m2 = Double_t(paras.GetParameterValue("m2"));
    double _m3 = Double_t(paras.GetParameterValue("m3"));
    //double locmax_sq = Double_t(paras.GetParameterValue("mb_"+name));
    //double locmin_sq = Double_t(paras.GetParameterValue("mb_"+name));
    unsigned int _subSysFlag = Double_t(paras.GetParameterValue("subSysFlag_"+name));
    double _inSpin = Double_t(paras.GetParameterValue("inSpin_"+name));
    double _outSpin1 = Double_t(paras.GetParameterValue("outSpin1_"+name));
    double _outSpin2 = Double_t(paras.GetParameterValue("outSpin2_"+name));

    double locmin_sq, locmax_sq, beta;

    switch(_subSysFlag){
      case 1:{ //reso in m23
        locmin_sq = s2min(_m23*_m23,_M,_m1,_m2,_m3);
        locmax_sq = s2max(_m23*_m23,_M,_m1,_m2,_m3);
        beta=acos((2.*_m13*_m13-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
        break;
      }
      case 2:{ //reso in m13
        locmin_sq = s1min(_m13*_m13,_M,_m1,_m2,_m3);
        locmax_sq = s1max(_m13*_m13,_M,_m1,_m2,_m3);
        beta=acos((2.*_m23*_m23-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
        break;
      }
      case 3:{ //reso in m12
        //return 1;
        locmin_sq = s1min(_m12*_m12,_M,_m1,_m3,_m2);
        locmax_sq = s1max(_m12*_m12,_M,_m1,_m3,_m2);
        beta=acos((2.*_m23*_m23-locmax_sq-locmin_sq)/(locmax_sq-locmin_sq));
        if(beta!=beta) return 1.;
        break;
      }
    }

    //double locmin_sq = s2min(_y*_y), locmax_sq = s2max(_y*_y);
    //if( _x*_x>locmax_sq || _x*_x<locmin_sq )
    //  return 0.;

    Spin j(_inSpin), m(_outSpin1), n(_outSpin2);
    return Wigner_d(j,m,n,beta);
  }

protected:
  std::string name;

};

class AmpSumTree {

public:
  /// Default Constructor (0x0)
  AmpSumTree(const double inM, const double inBr, const double in1
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
    m23(std::shared_ptr<DoubleParameter>(new DoubleParameter("m23", 1., m23_min, m23_max))), //ma
    m13(std::shared_ptr<DoubleParameter>(new DoubleParameter("m13", 1., m13_min, m13_max))), //mb
    m12(std::shared_ptr<DoubleParameter>(new DoubleParameter("m12", 1., m12_min, m12_max))), //mc
    x(std::shared_ptr<DoubleParameter>(new DoubleParameter("x", 0.))),
    totAmp("relBWsumAmplitude", "totAmp"){

    //AmpSumOfAmplitudes totAmp23("rbwAmp23", "totAmp");
    //------------Setup Tree---------------------
    myTree = std::shared_ptr<FunctionTree>(new FunctionTree());

    //----Strategies needed
    rbwStrat = std::shared_ptr<BreitWignerStrategy>(new BreitWignerStrategy());
    angdStrat = std::shared_ptr<WignerDStrategy>(new WignerDStrategy());
    multStrat = std::shared_ptr<MultAll>(new MultAll());
    addStrat = std::shared_ptr<AddAll>(new AddAll());

    //----Add Nodes
    myTree->createHead("Amplitude", addStrat); //A=Sum{Resos}

    //----Parameters needed
    //unsigned int numReso = ini.getResonances().size();
    for(std::vector<Resonance>::iterator reso=ini.getResonances().begin(); reso!=ini.getResonances().end(); reso++){
      Resonance tmp = (*reso);
      //setup RooVars
      mr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter(("m_{"+tmp.m_name+"}"), tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) );
      qr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter(("q_{"+tmp.m_name+"}"), tmp.m_breakup_mom) ) );
      gr.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter(("g_{"+tmp.m_name+"}"), tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) );
      std::complex<double> tmpintens = std::polar (tmp.m_strength, tmp.m_phase);
      rr.push_back( std::shared_ptr<ComplexParameter> (new ComplexParameter(("r_{"+tmp.m_name+"}"), tmpintens.real(),tmpintens.imag()) ) );
      //phir.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter(("phi_{"+tmp.m_name+"}"), tmp.m_phase) ) );

      aj.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter(("j_{"+tmp.m_name+"}"), tmp.m_spin) ) );
      am.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter(("m_{"+tmp.m_name+"}"), tmp.m_m) ) );
      an.push_back( std::shared_ptr<DoubleParameter> (new DoubleParameter(("n_{"+tmp.m_name+"}"), tmp.m_n) ) );

      //----Add Nodes
      unsigned int last = mr.size()-1;
      myTree->createNode("Reso_"+tmp.m_name, multStrat, "Amplitude"); //Reso=BW*c*AD
      myTree->createNode("RelBW_"+tmp.m_name, rbwStrat, "Reso_"+tmp.m_name); //BW
      myTree->createLeaf("Intens_"+tmp.m_name, rr[last], "Reso_"+tmp.m_name); //c
      myTree->createNode("AngD_"+tmp.m_name, angdStrat, "Reso_"+tmp.m_name); //AD
      //BW Par
      myTree->createLeaf("m0_"+tmp.m_name, mr[last], "RelBW_"+tmp.m_name); //m0
      myTree->createLeaf("x_"+tmp.m_name, x, "RelBW_"+tmp.m_name); //x
      myTree->createLeaf("ma", m23, "RelBW_"+tmp.m_name); //ma
      myTree->createLeaf("mb", m13, "RelBW_"+tmp.m_name); //mb
      myTree->createLeaf("spin_"+tmp.m_name, aj[last], "RelBW_"+tmp.m_name); //spin
      myTree->createLeaf("d_"+tmp.m_name, qr[last], "RelBW_"+tmp.m_name); //d
      myTree->createLeaf("resWidth_"+tmp.m_name, gr[last], "RelBW_"+tmp.m_name); //resWidth
      //AD Par
      myTree->createLeaf("m23", m23, "AngD_"+tmp.m_name); //ma
      myTree->createLeaf("m13", m13, "AngD_"+tmp.m_name); //mb
      myTree->createLeaf("m12", m12, "AngD_"+tmp.m_name); //mc
      myTree->createLeaf("M", M, "AngD_"+tmp.m_name); //M
      myTree->createLeaf("m1", m1, "AngD_"+tmp.m_name); //m1
      myTree->createLeaf("m2", m2, "AngD_"+tmp.m_name); //m2
      myTree->createLeaf("m3", m3, "AngD_"+tmp.m_name); //m3
     // unsigned int _subSysFlag = Double_t(paras.GetParameterValue("subSysFlag_"+name));
      //myTree->createLeaf("spin_"+tmp.m_name, aj[last], "AngD_"+tmp.m_name); //subSysFlag
      myTree->createLeaf("spin_"+tmp.m_name, aj[last], "AngD_"+tmp.m_name); //spin
      myTree->createLeaf("m_"+tmp.m_name, am[last], "AngD_"+tmp.m_name); //OutSpin 1
      myTree->createLeaf("n_"+tmp.m_name, an[last], "AngD_"+tmp.m_name); //OutSpin 2

      if(tmp.m_daugtherA==2 && tmp.m_daugtherB==3){
        std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
            tmp.m_name.c_str(), ma, *mr[last], *gr[last], *qr[last], 1, tmp.m_spin) );
        tmpbw->setDecayMasses(m2, m3);
        rbw.push_back(tmpbw);
        angd.push_back( std::shared_ptr<AmpWigner> (new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
            mc, ma, mb, 1, *aj[last], *am[last], *an[last]) ) );
        totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last), angd.at(last));
      }else if(tmp.m_daugtherA==1 && tmp.m_daugtherB==3){
        std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
            tmp.m_name.c_str(), mb, *mr[last], *gr[last], *qr[last], 2, tmp.m_spin) );
        tmpbw->setDecayMasses(m1, m3);
        rbw.push_back(tmpbw);
        angd.push_back( std::shared_ptr<AmpWigner> (new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
            mc, ma, mb, 2, *aj[last], *am[last], *an[last]) ) );
        totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last), angd.at(last));
      }else if(tmp.m_daugtherA==1 && tmp.m_daugtherB==2){
        std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(tmp.m_name.c_str(),
            tmp.m_name.c_str(), mc, *mr[last], *gr[last], *qr[last], 3, tmp.m_spin) );
        tmpbw->setDecayMasses(m1, m2);
        rbw.push_back(tmpbw);
        angd.push_back( std::shared_ptr<AmpWigner> (new AmpWigner(("a_{"+tmp.m_name+"}").c_str(), ("a_{"+tmp.m_name+"}").c_str(),
            mc, ma, mb, 3, *aj[last], *am[last], *an[last]) ) );
        totAmp.addBW(rbw.at(last), rr.at(last), phir.at(last), angd.at(last));
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

  virtual const double volume(){
    return -1.*(m23_sq_min-m23_sq_max)*(m13_sq_min-m13_sq_max)*(m12_sq_min-m12_sq_max);
  }

  virtual const double intensity(std::vector<double>& x, ParameterList& par){
    //TODO: check x exception
    if(x.size()!=3) return 0;

    ma.setVal(x[0]); mb.setVal(x[1]); mc.setVal(x[2]);

    if( par.GetNDouble()>mr.size() ){
        std::cout << "Error: Parameterlist doesn't match model!! " << par.GetNDouble() << std::endl; //TODO: exception
      return 0;
    }

    for(unsigned int i=0; i<mr.size(); i++){
      mr[i]->setVal(par.GetDoubleParameter(i)->GetValue());
      //phir[i]->setVal(par.GetDoubleParameter(2*i+1).GetValue());
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

    //add strength and phases of the used amplitudes
    for(unsigned int i=0; i<rr.size();i++){


      if(rr[i]->hasError()) //TODO: check bounds
        outPar.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter(mr[i]->GetName(),mr[i]->getVal(), mr[i]->getMin(), mr[i]->getMax(), mr[i]->getError())));
      else
        outPar.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter(mr[i]->GetName(),mr[i]->getVal(), 0.1)));

    }

    return true;
  }

  virtual std::shared_ptr<FunctionTree> functionTree(ParameterList& outPar) {
    if(outPar.GetNParameter()>0) return NULL;
    fillStartParVec(outPar);

    return myTree;
  }

  /** Destructor */
  virtual ~AmpSumTree(){

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

  std::shared_ptr<DoubleParameter> ma;
  std::shared_ptr<DoubleParameter> mb;
  std::shared_ptr<DoubleParameter> mc;
  std::shared_ptr<DoubleParameter> x;

  //AmpSumOfAmplitudes totAmp;
  std::shared_ptr<FunctionTree> myTree;

  //Resonance Variables
  std::vector<std::shared_ptr<DoubleParameter> > mr;
  std::vector<std::shared_ptr<DoubleParameter> > qr;
  std::vector<std::shared_ptr<DoubleParameter> > gr;
  std::vector<std::shared_ptr<ComplexParameter> > rr;
  //std::vector<std::shared_ptr<DoubleParameter> > phir;
  std::vector<std::shared_ptr<DoubleParameter> > aj;
  std::vector<std::shared_ptr<DoubleParameter> > am;
  std::vector<std::shared_ptr<DoubleParameter> > an;

  std::shared_ptr<BreitWignerStrategy> rbwStrat;
  std::shared_ptr<WignerDStrategy> angdStrat;
  std::shared_ptr<MultAll> multStrat;
  std::shared_ptr<AddAll> addStrat;

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
