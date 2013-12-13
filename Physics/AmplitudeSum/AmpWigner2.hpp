//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel
//     Peter Weidenkaff
//-------------------------------------------------------------------------------

//! Angular distribution based on WignerD functions
/*!
 * @file AmpWigner2.hpp
 *\class AmpWigner2
 *The helicity angle for sub system \_subSys is calculated and the value of the WignerD function is returned
 */

#ifndef AMPWIGNER2
#define AMPWIGNER2

#include <vector>

#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Core/DataPoint.hpp"

using namespace std;

class AmpWigner2{
public:
  AmpWigner2( unsigned int subSys, unsigned int resSpin);

  AmpWigner2(const AmpWigner2&, const char*);

  virtual ~AmpWigner2() {};

  virtual double evaluate(dataPoint& point) const;
  void setDecayMasses(double, double, double, double);

protected:

  virtual void initialise();

  unsigned int _resSpin;
  unsigned int _subSys;

  double _M;
  double _m1;
  double _m2;
  double _m3;

 unsigned int _spinM;
 unsigned int _spin1;
 unsigned int _spin2;
 unsigned int _spin3;
private:

};

#endif
