//****************************************************************************
// Abstract base class for dynamical functions.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Abstract base class for dynamical functions.

#ifndef AMP_ABS_DYNAMICAL_FUNCTION
#define AMP_ABS_DYNAMICAL_FUNCTION

#include "TObject.h"
#include "RooComplex.h"
#include "RooAbsArg.h"

#include <vector>
#include <complex>

class AmpAbsDynamicalFunction : public RooAbsArg  {
public:
  AmpAbsDynamicalFunction(const char *name, const char *title);

  AmpAbsDynamicalFunction(const AmpAbsDynamicalFunction&, const char*);

  virtual ~AmpAbsDynamicalFunction();


  virtual void initialise() = 0; 
  virtual RooComplex evaluate() const = 0;
//  virtual double evaluate(double x[],int dim, void * param) const = 0;//used for MC integration
  virtual double getSpin() = 0;
  virtual bool isSubSys(const unsigned int) const = 0;
  virtual double integral() const = 0;
 

  // the following are needed by the RooAbsArg interface, but not yet 
  // implemented

  virtual TObject*  clone (const char *newname) const = 0 ;  


protected:
//  virtual double evaluateAngle() const = 0;

private:

  //ClassDef(AmpAbsDynamicalFunction,1) // Abstract dynamical function

};

#endif
