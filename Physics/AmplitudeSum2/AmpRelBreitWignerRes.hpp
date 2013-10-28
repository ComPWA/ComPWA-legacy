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

//#include "TObject.h"
//#include "TString.h"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpKinematics.hpp"
#include "Physics/AmplitudeSum/AmpWigner.hpp"

class AmpRelBreitWignerRes : public AmpAbsDynamicalFunction, public AmpKinematics {
public:

  AmpRelBreitWignerRes(const char *name,
		       DoubleParameter& _resMass, DoubleParameter& _resWidth,
		       double& _radius,
		       int _subsys,
               int resSpin, int m, int n
               ) ;
  AmpRelBreitWignerRes(const AmpRelBreitWignerRes&, const char*);
  AmpRelBreitWignerRes(const AmpRelBreitWignerRes&);

  virtual ~AmpRelBreitWignerRes();

  virtual void initialise();
  virtual std::complex<double> evaluate() const ;
  void setDecayMasses(double, double, double, double);
//  double integral(ParameterList& par) const;
//  double getMaximum() const{return 1;};
  double integral() const;
  double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
  inline virtual bool isSubSys(const unsigned int subSys) const{ return (subSys==_subSys); };

protected:
  virtual double eval(double x[],int dim, void * param) const;//used for MC integration

  DoubleParameter _resWidth;
  AmpWigner _wignerD;

private:
  //ClassDef(AmpRelBreitWignerRes,1) // Relativistic Breit-Wigner resonance model

};

#endif
