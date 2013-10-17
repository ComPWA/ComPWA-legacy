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
#include "Physics/AmplitudeSum/AmpKinematics.hpp"
#include "Physics/AmplitudeSum/AmpWigner.hpp"

class AmpRelBreitWignerRes : public AmpAbsDynamicalFunction, public AmpKinematics {
public:

  AmpRelBreitWignerRes(const char *name, const char *title,
		       RooAbsReal& _resMass, RooAbsReal& _resWidth,
		       RooAbsReal& _radius,
		       int _subsys,
               Int_t resSpin,
               Int_t m,
               Int_t n
               ) ;
  AmpRelBreitWignerRes(const char *name, const char *title,
		       RooAbsReal& _x13, ///  second DP variable
		       RooAbsReal& _x23, ///  mass at which to evaluate RBW
		       RooAbsReal& _resMass, RooAbsReal& _resWidth,
		       RooAbsReal& _q0,
		       int _subsys,
               Int_t resSpin,
               Int_t m,
               Int_t n
               ) ;

  AmpRelBreitWignerRes(const AmpRelBreitWignerRes&, const char*);
  AmpRelBreitWignerRes(const AmpRelBreitWignerRes&);

  virtual ~AmpRelBreitWignerRes();

  virtual void initialise();
  virtual RooComplex evaluate() const ;
//  virtual double evaluate(double x[],int dim, void * param) const;//used for MC integration
  void setDecayMasses(double, double, double, double);
//  double integral(ParameterList& par) const;
//  double getMaximum() const{return 1;};
  double integral() const {return 1;};
	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction

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
  RooRealProxy _x13;
  RooRealProxy _x23;
  RooRealProxy _resWidth;
//  unsigned int _subSys;
  AmpWigner _wignerD;

//  virtual double evaluateAngle() const {};
private:

  //ClassDef(AmpRelBreitWignerRes,1) // Relativistic Breit-Wigner resonance model

};

#endif
