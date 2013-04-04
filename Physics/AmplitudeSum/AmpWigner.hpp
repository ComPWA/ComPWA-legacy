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
		       RooAbsReal& _m12, RooAbsReal& _m23, RooAbsReal& _m13,
                       const int _subSysFlag,
		       RooAbsReal& _inSpin, RooAbsReal& _outSpin1,
		       RooAbsReal& _outSpin2) ; 

  AmpWigner(const AmpWigner&, const char*);

  virtual ~AmpWigner();

  virtual inline bool hasDist(){return toEvaluate;};
  virtual void initialise();
  virtual double evaluate() const;

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

protected:
  RooRealProxy _m12; //m12
  RooRealProxy _m23; //m23
  RooRealProxy _m13; //m13

  int _subSysFlag;

  RooRealProxy _inSpin; 
  RooRealProxy _outSpin1;
  RooRealProxy _outSpin2;

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
