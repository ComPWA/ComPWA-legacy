//****************************************************************************
// Class for defining Wigner_d angular distributions
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining Wigner_d angular distributions

#ifndef AMPWIGNER
#define AMPWIGNER

#include <vector>

#include "RooRealProxy.h"
#include "RooAbsReal.h"

using namespace std;

class AmpWigner : public RooAbsArg {
public:

  AmpWigner();

  AmpWigner(const char *name, const char *title,
		       RooAbsReal& _m13, RooAbsReal& _m23,
		       UInt_t spin, UInt_t m, UInt_t n) ;

  AmpWigner(const AmpWigner&, const char*);

  virtual ~AmpWigner();

  virtual inline bool hasDist(){return toEvaluate;};
  virtual void initialise();
  virtual double evaluate() const;
  void setDecayMasses(double, double, double, double);

  // the following are needed by the RooAbsArg interface, but not yet 
  // implemented

  virtual TObject*  clone (const char *newname) const ;
  virtual Bool_t readFromStream(std::istream&, Bool_t, Bool_t);
  virtual void writeToStream(std::ostream&, Bool_t) const;
  virtual Bool_t operator==(const RooAbsArg&);
  virtual void syncCache(const RooArgSet*);
  virtual void copyCache(const RooAbsArg*, Bool_t, Bool_t);
  virtual void attachToTree(TTree&, Int_t);
  virtual void attachToVStore(RooVectorDataStore&);
  virtual void setTreeBranchStatus(TTree&, Bool_t);
  virtual void fillTreeBranch(TTree&);
  virtual RooAbsArg* createFundamental(const char*) const ;
  virtual Bool_t isIdentical(const RooAbsArg&, Bool_t);

protected:
  RooRealProxy _m23; //m23
  RooRealProxy _m13; //m13


  UInt_t _inSpin;
  UInt_t _outSpin1;
  UInt_t _outSpin2;

  Double_t _M;
  Double_t _m1;
  Double_t _m2;
  Double_t _m3;

  bool toEvaluate;

  double lambda(double x, double y, double z)const;
  Double_t s1max(Double_t, Double_t, Double_t, Double_t, Double_t)const;
  Double_t s1min(Double_t, Double_t, Double_t, Double_t, Double_t)const;
  Double_t s2max(Double_t, Double_t, Double_t, Double_t, Double_t)const;
  Double_t s2min(Double_t, Double_t, Double_t, Double_t, Double_t)const;
  Double_t s3max(Double_t, Double_t, Double_t, Double_t, Double_t)const;
  Double_t s3min(Double_t, Double_t, Double_t, Double_t, Double_t)const;

private:

  //ClassDef(AmpWigner,1) // Wigner_d angular distribution

};

#endif
