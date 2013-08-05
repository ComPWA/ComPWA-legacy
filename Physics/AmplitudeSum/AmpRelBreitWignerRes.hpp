//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef AMP_REL_BREIT_WIGNER_RES
#define AMP_REL_BREIT_WIGNER_RES

#include <vector>

#include "TObject.h"
#include "TString.h"
#include "RooComplex.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooRealProxy.h"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"

class AmpRelBreitWignerRes : public AmpAbsDynamicalFunction  {
public:

  AmpRelBreitWignerRes(const char *name, const char *title,
		       RooAbsReal& _x, ///  mass at which to evaluate RBW
		       RooAbsReal& _resMass, RooAbsReal& _resWidth,
		       RooAbsReal& _q0,
		       int _subsys,
               Int_t resSpin) ;

  AmpRelBreitWignerRes(const AmpRelBreitWignerRes&, const char*);
  AmpRelBreitWignerRes(const AmpRelBreitWignerRes&);


  virtual ~AmpRelBreitWignerRes();

  void setDecayMasses(double, double);
  double getSpin(){return _spin;};

  double q0() const;
  double q()  const;
  double BLprime2() const;
  double F(double) const;
  
  virtual void initialise();
  virtual RooComplex evaluate() const ;

  // the following are needed by the RooAbsArg interface, but not yet 
  // implemented

  virtual TObject*  clone (const char *newname) const ;
  virtual Bool_t readFromStream(std::istream&, Bool_t, Bool_t);
  virtual void writeToStream(std::ostream&, Bool_t) const;
  virtual Bool_t operator==(const RooAbsArg&);
  virtual void syncCache(const RooArgSet*);
  virtual void copyCache(const RooAbsArg*, Bool_t, Bool_t);
  virtual void attachToTree(TTree&, Int_t);
  virtual void setTreeBranchStatus(TTree&, Bool_t);
  virtual void attachToVStore(RooVectorDataStore&);
  virtual void fillTreeBranch(TTree&);
  virtual RooAbsArg* createFundamental(const char*) const ;
  virtual bool isIdentical(const RooAbsArg&, Bool_t);

  inline virtual bool isSubSys(const unsigned int subSys)const{return (subSys==_subSys);};

protected:
  RooRealProxy _x;

  RooRealProxy _m0;
  RooRealProxy _resWidth;
  RooRealProxy _d;
  unsigned int _subSys;
  int _spin;

  // masses of decay particles for this resonance
  double _ma;
  double _mb;


private:

  //ClassDef(AmpRelBreitWignerRes,1) // Relativistic Breit-Wigner resonance model

};

#endif
