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
#include "RooAbsArg.h"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"

AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(const char *name, const char *title) :
  RooAbsArg (name, title) 
{

}

AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(const AmpAbsDynamicalFunction& other, const char* newname) :
  RooAbsArg(other, newname) 
{
}

AmpAbsDynamicalFunction::~AmpAbsDynamicalFunction() 
{
}


TObject* AmpAbsDynamicalFunction::clone(const char *newname)  const
{
  return 0;
}
