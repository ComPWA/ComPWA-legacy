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

//#include "TObject.h"
//#include "RooComplex.h"
//#include "RooAbsArg.h"
#include "Core/Parameter.hpp"
#include <vector>
#include <complex>

class AmpAbsDynamicalFunction{// : public Strategy{
public:
  AmpAbsDynamicalFunction(const char *name);

  AmpAbsDynamicalFunction(const AmpAbsDynamicalFunction&, const char*);

  virtual ~AmpAbsDynamicalFunction();

  //! Implementation of interface for streaming info about the strategy
  virtual const std::string to_str() const {
	  return (_name+"rel Breit-Wigner with Wigner_D");
  }

  //! Implementation of interface for executing a strategy
 /* virtual std::shared_ptr<AbsParameter> execute(const ParameterList& paras){
	  std::complex<double> result;

    result = evaluateTree(paras);

    //ParameterList out;
    // out.AddParameter(DoubleParameter("AddAllResult",result));
    return std::shared_ptr<AbsParameter>(new ComplexParameter("AddAllResult",result));
  };*/


  virtual void initialise() = 0; 
  virtual std::complex<double> evaluate() const = 0;
  virtual std::complex<double> evaluateAmp() const = 0;
  /*virtual std::complex<double> evaluateTree(const ParameterList& paras) const {
	  //TODO: pure virtual
	  std::complex<double> blanko;
	  return blanko;
  }*/
  virtual double evaluateWignerD() const = 0;

  virtual double getSpin() = 0;
  virtual bool isSubSys(const unsigned int) const = 0;
  virtual double evaluate(double x[],size_t dim) const;//used for MC integration
	//! Calculation integral |resonance|^2
  virtual double integral() const;

  virtual std::string GetName(){ return _name; };
  virtual std::string GetTitle(){ return GetName(); };
  virtual double GetNormalization(){ return _norm; };
  virtual void SetNormalization(double n){ _norm=n; };
 
protected:
  std::string _name;
  double _norm;

private:
  //ClassDef(AmpAbsDynamicalFunction,1) // Abstract dynamical function

};

#endif
