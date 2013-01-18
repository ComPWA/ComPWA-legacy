#include "RooAbsArg.h"
#include "AmpAbsDynamicalFunction.hpp"

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
