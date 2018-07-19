/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtDalitzReso.hh,v 1.1 2009-03-16 16:50:49 robbep Exp $
 *
 * Description:
 *   Class to compute Dalitz amplitudes based on many models that cannot be
 *     handled with EvtResonance.
 *
 * Modification history:
 *   Jordi Garra Tic�     2008/07/03         File created
 *****************************************************************************/

#ifndef __EVTDALITZRESO_HH__
#define __EVTDALITZRESO_HH__

#include <string>
#include <vector>
#include <map>

#include "EvtComplex.hh"
#include "EvtCyclic3.hh"
#include "EvtSpinType.hh"
#include "EvtTwoBodyVertex.hh"
#include "EvtDalitzPoint.hh"
// #include "EvtDecayAmp.hh"
#include "EvtBlattWeisskopf.hh"
#include "EvtFlatte.hh"

using std::vector;
using std::map;

class EvtBlattWeisskopf;

class EvtDalitzReso
{
public:
  // Numerator type
  enum NumType { NBW            = 0 , RBW_ZEMACH        = 1 , RBW_KUEHN  = 2 , RBW_CLEO        = 3 ,
		 RBW_ZEMACH2    = 4 , GS_CLEO           = 5 , K_MATRIX   = 6 , RBW_CLEO_ZEMACH = 7 ,
		 GS_CLEO_ZEMACH = 8 , LASS              = 9 , K_MATRIX_I = 10, K_MATRIX_II     = 11,
		 GAUSS_CLEO     = 12, GAUSS_CLEO_ZEMACH = 13, FLATTE = 14, NON_RES = 15 };

  // Coupling type
  //  ChgPion : pi+ pi-
  //  NeuPion : pi0 pi0
  //  Pion    : 0.5*[(pi+ pi-) + (pi0 pi0)]
  //  ChgKaon : K+ K-
  //  NeuKaon : K0 K0
  //  Kaon    : 0.5*[(K+ K-) + (K0 K0)]
  //  EtaPion : eta pi0
  enum CouplingType {Undefined=0,PicPic=1,PizPiz,PiPi,KcKc,KzKz,KK,EtaPic,EtaPiz,PicPicKK,WA76};

  EvtDalitzReso() : _typeN(NON_RES) {};

  EvtDalitzReso(std::string name) : _Name(name), _typeN(NON_RES) {};
  
  EvtDalitzReso(std::string name, const EvtDalitzPlot& dp, EvtCyclic3::Pair pairAng, EvtCyclic3::Pair pairRes, 
		EvtSpinType::spintype spin, double m0, double g0, NumType typeN, double f_b=0.0, double f_d=1.5);

  EvtDalitzReso(std::string name, const EvtDalitzPlot& dp, EvtCyclic3::Pair pairAng, EvtCyclic3::Pair pairRes, 
		EvtSpinType::spintype spin, double m0, double g0, NumType typeN,
		double m0_mix, double g0_mix, double delta_mix, EvtComplex amp_mix);

  EvtDalitzReso(std::string name, const EvtDalitzPlot& dp, EvtCyclic3::Pair pairAng, EvtCyclic3::Pair pairRes, 
		EvtSpinType::spintype spin, double m0, NumType typeN, double g1, double g2, CouplingType coupling2);

  // K-matrix
  EvtDalitzReso(std::string name, const EvtDalitzPlot& dp, EvtCyclic3::Pair pairRes, std::string nameIndex, NumType typeN,
		EvtComplex fr12prod, EvtComplex fr13prod, EvtComplex fr14prod, EvtComplex fr15prod, double s0prod);

  // K-matrix FIT
  EvtDalitzReso(std::string name, const EvtDalitzPlot& dp, NumType typeN); 
  
  // LASS
  EvtDalitzReso(std::string name, const EvtDalitzPlot& dp, EvtCyclic3::Pair pairRes, double m0, double g0,
		double a, double r, double B, double phiB, double R, double phiR);

  //Flatte
  EvtDalitzReso(std::string name, const EvtDalitzPlot& dp, EvtCyclic3::Pair pairRes, double m0);

  EvtDalitzReso(const EvtDalitzReso& other);

  ~EvtDalitzReso();

  EvtDalitzReso* clone() const { return new EvtDalitzReso(*this); }

  EvtComplex evaluate(const EvtDalitzPoint& p);
  vector< EvtComplex > evaluateKMat(double s);
  
  vector< EvtComplex > evaluateTOneOne(double s);
  
  // Setter Functions
  void set_fd( double R ) { _vd.set_f( R ); }
  void set_fb( double R ) { _vb.set_f( R ); }

  void set_Mass( double Mass ) { _m0=Mass; }
  void set_Gamma( double Gamma ) { _g0=Gamma; }

  double get_Mass() const { return _m0; }
  double get_Gamma() const { return _g0; }
  
  void set_a( double a );
  void set_r( double r ); 
  void set_B( double B ); 
  void set_phiB( double phiB );
  void set_R( double Rneu );
  void set_phiR( double phiR );
  
  // Getter Functions
  std::string get_Name() const { return _Name; }
  NumType get_Type() const { return _typeN; }
  
  double get_a() const { return _a; }
  double get_r() const {  return _r; }
  double get_B() const {  return _Blass;  } 
  double get_phiB() const {  return _phiB; }
  double get_R() const {  return _R; }
  double get_phiR() const {  return _phiR; }
  
//   void Print_Mass( void ) { std::cout<<"Mass: "<<_m0<<std::endl; }
//   void Print_Gamma( void ) { std::cout<<"Gamma: "<<_g0<<std::endl; }
  
//   void Print_All( void ) {
//   std::cout<<"-----------------------"		<<std::endl;  
//   std::cout<<"Mass: "	<<_m0	       		<<std::endl;
//   std::cout<<"Gamma: "<<_g0           		<<std::endl;  
//   std::cout<<"pairAng: "<<c_str(_pairAng)      <<std::endl; 
//   std::cout<<"pairRes: "<<c_str(_pairRes)	<<std::endl;
//   std::cout<<"Spin: "<<EvtSpinType::getSpin2(_spin)		<<std::endl; 
//   std::cout<<"Type: "<<_typeN			<<std::endl;   
//   }
// 
//   void Print_KMat( void ) {
//   std::cout<<"-----------------------"<<std::endl;   
//   std::cout<<"_pairRes: "<<c_str(_pairRes) <<std::endl; 
//   std::cout<<"_fr12prod: "<<_fr12prod<<std::endl;
//   std::cout<<"_fr13prod: "<<_fr13prod<<std::endl;
//   std::cout<<"_fr14prod: "<<_fr14prod<<std::endl;
//   std::cout<<"_fr15prod: "<<_fr15prod<<std::endl;
//   std::cout<<"_s0prod: "<<_s0prod<<std::endl; 
//   std::cout<<"_kmatrix_index: "<<_kmatrix_index<<std::endl;
//   std::cout<<"Type: "<<_typeN			<<std::endl;   
//   }
    
  void addFlatteParam(const EvtFlatteParam& param) { _flatteParams.push_back(param); }

  EvtComplex Fvector( double s, int index );

  EvtComplex CalcdeltaR(const EvtDalitzPoint& x);
  EvtComplex evaluateLASS(const EvtDalitzPoint& x,const EvtComplex& deltaR_comp);  
  
  EvtComplex lass(double s);
  
  EvtComplex CalcGounarisSakurai(double m);
  
private:
  EvtComplex psFactor(double& ma, double& mb, double& m);
  EvtComplex psFactor(double& ma1, double& mb1, double& ma2, double& mb2, double& m);
  EvtComplex propGauss(const double& m0, const double& s0, const double& m);
  EvtComplex propBreitWigner(const double& m0, const double& g0, const double& m); 
  EvtComplex propBreitWignerRel(const double& m0, const double& g0, const double& m);
  EvtComplex propBreitWignerRel(const double& m0, const EvtComplex& g0, const double& m);
  EvtComplex propBreitWignerRelCoupled(const double& m0, const EvtComplex& g1, const EvtComplex& g2, const double& m); 
  EvtComplex propGounarisSakurai(const double& m0, const double& g0, const double& k0, 
	                         const double& m, const double& g, const double& k);
  inline double GS_f(const double& m0, const double& g0, const double& k0, const double& m, const double& k);
  inline double GS_h(const double& m, const double& k);
  inline double GS_dhods(const double& m0, const double& k0); 
  inline double GS_d(const double& m0, const double& k0); 

  EvtComplex numerator(const EvtDalitzPoint& p, const EvtTwoBodyKine& vb, const EvtTwoBodyKine& vd);
  double angDep(const EvtDalitzPoint& p);
  EvtComplex mixFactor(EvtComplex prop, EvtComplex prop_mix);
  
  EvtComplex flatte(const double& m);

  inline EvtComplex sqrtCplx(double in) { return (in > 0) ? EvtComplex(sqrt(in), 0) : EvtComplex(0, sqrt(-in)); }

    // ResoName
  std::string _Name;
  
  // Dalitz plot
  EvtDalitzPlot _dp; 

  // Pairing indices:
  EvtCyclic3::Pair _pairAng;    // angular  
  EvtCyclic3::Pair _pairRes;    // resonance

  // Spin
  EvtSpinType::spintype _spin;                                  

  // Numerator type
  NumType _typeN;

  // Nominal mass and width
  double _m0,_g0; 

  // Vertices
  EvtTwoBodyVertex _vb;
  EvtTwoBodyVertex _vd;

  // Daughter masses
  double _massFirst,_massSecond;

  // variables for electromagnetic mass mixing 
  double _m0_mix,_g0_mix,_delta_mix;
  EvtComplex _amp_mix;   

  // variables for coupled Breit-Wigner
  double _g1,_g2;
  CouplingType _coupling2;

  // variables for Blatt-Weisskopf form factors
  double _f_b, _f_d;

  // K-matrix 
  int _kmatrix_index;
  EvtComplex _fr11prod, _fr12prod,_fr13prod,_fr14prod,_fr15prod;
  double _s0prod;

  // LASS
  double _a;
  double _r;
  double _Blass;
  double _phiB;
  double _R;
  double _phiR;

  // Flatte
  std::vector<EvtFlatteParam> _flatteParams;
  
//   EvtComplex beta1, beta2, beta3, beta4, beta5;
  
};

#endif

