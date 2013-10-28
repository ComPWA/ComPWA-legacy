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
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpGausRes.hpp"

AmpGausRes::AmpGausRes(const char *name,
					   DoubleParameter& resMass, DoubleParameter& resWidth,
					   int subSys) :
  AmpAbsDynamicalFunction(name),
  _mR(resMass), _resWidth(resWidth),
  _subSys(subSys)
{
  initialise();
}


AmpGausRes::AmpGausRes(const AmpGausRes& other, const char* newname) :
  AmpAbsDynamicalFunction(other, newname),
  _mR(other._mR),
  _resWidth(other._resWidth),
  _subSys(other._subSys)
{
  initialise();
}

AmpGausRes::AmpGausRes(const AmpGausRes& other) :
  AmpAbsDynamicalFunction(other),
  _mR(other._mR),
  _resWidth(other._resWidth),
  _subSys(other._subSys)
{
  initialise();
}

AmpGausRes::~AmpGausRes() 
{
}

void AmpGausRes::initialise() 
{
}   

std::complex<double> AmpGausRes::evaluate() const {


  double m0 = _mR.GetValue();
  double width = _resWidth.GetValue();
//  double m  = Double_t(_x);
	double m = dataPoint::instance()->getM(_subSys);
  
  std::complex<double> gaus (exp(-1*(m-m0)*(m-m0)/width/width/2.),0);

  return gaus;
}
