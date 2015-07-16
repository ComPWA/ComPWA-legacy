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

#include <vector>
#include <complex>

#include "Core/Parameter.hpp"
#include "Core/DataPoint.hpp"

class AmpAbsDynamicalFunction
{

public:
  AmpAbsDynamicalFunction(const char *name, int nCalls=30000);

  virtual ~AmpAbsDynamicalFunction();

  //! Implementation of interface for streaming info about the strategy
  virtual const std::string to_str() const { return (_name); }
  //! value of resonance at \param point
  virtual std::complex<double> evaluate(dataPoint& point) = 0;
  //! value of dynamical amplitude at \param point
  virtual std::complex<double> evaluateAmp(dataPoint& point) = 0;

  //! Get resonance spin
  virtual double getSpin() = 0;
  //! Calculation integral |dynamical amplitude|^2
  virtual double integral();
  //! Calculation integral |dynamical amplitude * WignerD|^2
  virtual double totalIntegral() const;

  virtual std::string GetName(){ return _name; };
  virtual std::string GetTitle(){ return GetName(); };
  //! Get current normalization
  virtual double GetNormalization();
  //! Set normalization manually. Setting to values <0 disables normalization
  virtual void SetNormalization(double n){ _norm=n; };
 
  virtual void SetModified() { modified=1; }

protected:
  std::string _name;
  int _nCalls;

private:
  double _norm;
  bool modified;

};

#endif
