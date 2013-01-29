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

class AmpAbsDynamicalFunction : public RooAbsArg  {
public:
  AmpAbsDynamicalFunction(const char *name, const char *title);

  AmpAbsDynamicalFunction(const AmpAbsDynamicalFunction&, const char*);

  virtual ~AmpAbsDynamicalFunction();

  virtual void initialise() = 0; 
  virtual RooComplex evaluate() = 0;
  virtual bool isSubSys(const unsigned int) = 0;
 

  // the following are needed by the RooAbsArg interface, but not yet 
  // implemented

  virtual TObject*  clone (const char *newname) const = 0 ;  


protected:

private:

  //ClassDef(AmpAbsDynamicalFunction,1) // Abstract dynamical function

};

#endif
