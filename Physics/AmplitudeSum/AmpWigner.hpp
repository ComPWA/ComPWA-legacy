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
#include "Physics/DPKinematics/DPKinematics.hpp"
#include "Physics/DPKinematics/DataPoint.hpp"

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
  virtual double evaluate() const;
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

private:
  //ClassDef(AmpWigner,1) // Wigner_d angular distribution

};

#endif
