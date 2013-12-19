//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//     Peter Weidenkaff -  assignment of final state particle masses
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining Wigner_d angular distributions
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining Wigner_d angular distributions

#ifndef AMPWIGNER
#define AMPWIGNER

#include <vector>

//#include "RooRealProxy.h"
//#include "RooAbsReal.h"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Core/DataPoint.hpp"

using namespace std;

class AmpWigner{
public:

  AmpWigner();
  AmpWigner(const char *name,
		       unsigned int spin, unsigned int m, unsigned int n, unsigned int subSys) ;

  AmpWigner(const AmpWigner&, const char*);

  virtual ~AmpWigner();

  virtual inline bool hasDist(){return toEvaluate;};
  virtual void initialise();
  virtual double evaluate(dataPoint& point) const;
  void setDecayMasses(double, double, double, double);

protected:
  unsigned int _inSpin;
  unsigned int _outSpin1;
  unsigned int _outSpin2;
  unsigned int _subSys;

  double _M;
  double _m1;
  double _m2;
  double _m3;

  bool toEvaluate;

  double lambda(double x, double y, double z)const;
  double s1max(double, double, double, double, double)const;
  double s1min(double, double, double, double, double)const;
  double s2max(double, double, double, double, double)const;
  double s2min(double, double, double, double, double)const;
  double s3max(double, double, double, double, double)const;
  double s3min(double, double, double, double, double)const;

  double cosTheta(double, double, double)const;
  double qin(double)const;
  double qout(double)const;

private:
  //ClassDef(AmpWigner,1) // Wigner_d angular distribution

};

#endif
