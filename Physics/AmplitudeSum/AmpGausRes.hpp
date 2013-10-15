//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef AMP_GAUS_RES
#define AMP_GAUS_RES

#include <vector>

#include "TObject.h"
#include "TString.h"
#include "RooComplex.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooRealProxy.h"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"

using namespace std;

class AmpGausRes : public AmpAbsDynamicalFunction  {
public:

  AmpGausRes(const char *name, const char *title,
		       RooAbsReal& _x, ///  mass at which to evaluate RBW
		       RooAbsReal& _resMass, RooAbsReal& _resWidth,
		       int _subsys) ; 

  AmpGausRes(const AmpGausRes&, const char*);
  AmpGausRes(const AmpGausRes&);

  ~AmpGausRes();

  virtual void initialise();
  virtual RooComplex evaluate()const;
  virtual double evaluateAngle() const {};

  double getSpin(){return 0;};
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

  inline virtual bool isSubSys(const unsigned int subSys)const{return (subSys==_subSys);};

protected:
  RooRealProxy _x;

  RooRealProxy _m0;
  RooRealProxy _resWidth;

  unsigned int _subSys;


private:

  //ClassDef(AmpGausRes,1) // Relativistic Breit-Wigner resonance model

};

#endif
