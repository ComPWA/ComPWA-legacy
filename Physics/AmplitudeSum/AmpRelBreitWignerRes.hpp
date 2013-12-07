//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#ifndef AMP_REL_BREIT_WIGNER_RES
#define AMP_REL_BREIT_WIGNER_RES

#include <vector>

//#include "TObject.h"
//#include "TString.h"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpKinematics.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Physics/AmplitudeSum/AmpWigner.hpp"

class AmpRelBreitWignerRes : public AmpAbsDynamicalFunction, public AmpKinematics {
public:

	AmpRelBreitWignerRes(const char *name,
			DoubleParameter& _resMass, DoubleParameter& _resWidth,
			double& _radius,
			int _subsys,
			int resSpin, int m, int n
	) ;
	AmpRelBreitWignerRes(const AmpRelBreitWignerRes&, const char*);
	AmpRelBreitWignerRes(const AmpRelBreitWignerRes&);

	virtual ~AmpRelBreitWignerRes();

	//  double operator() (double *x, size_t dim, void*);

	virtual void initialise();
	virtual std::complex<double> evaluate() const ;
	virtual std::complex<double> evaluateAmp() const;
	virtual double evaluateWignerD() const { return _wignerD.evaluate(); };
	//  virtual double eval(double x[],size_t dim, void *param) const;//used for MC integration
	//  double (*eval2)(double x[],size_t dim, void *param);//used for MC integration

	void setDecayMasses(double, double, double, double);
	//  double getMaximum() const{return 1;};
	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	inline virtual bool isSubSys(const unsigned int subSys) const{ return (subSys==_subSys); };

protected:

	DoubleParameter _resWidth;
	AmpWigner2 _wignerD;
//	AmpWigner _wignerD;

private:

};

#endif
