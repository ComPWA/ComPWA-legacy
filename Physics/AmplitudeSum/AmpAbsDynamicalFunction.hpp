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
  virtual RooComplex evaluate() const = 0;
  virtual bool isSubSys(const unsigned int) const = 0;
 

  // the following are needed by the RooAbsArg interface, but not yet 
  // implemented

  virtual TObject*  clone (const char *newname) const = 0 ;  


protected:

private:

  //ClassDef(AmpAbsDynamicalFunction,1) // Abstract dynamical function

};

#endif
