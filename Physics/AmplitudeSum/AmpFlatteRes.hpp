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

#ifndef AMP_FLATTE_RES
#define AMP_FLATTE_RES

#include <vector>

//#include "TObject.h"
//#include "TString.h"
//#include "RooComplex.h"
//#include "RooAbsReal.h"
//#include "RooAbsArg.h"
//#include "RooRealProxy.h"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpKinematics.hpp"
#include "Physics/AmplitudeSum/AmpWigner.hpp"

using namespace std;

class AmpFlatteRes : public AmpAbsDynamicalFunction, public AmpKinematics {
public:
	AmpFlatteRes(const char *name,
			DoubleParameter& resMass, DoubleParameter& resWidth,
			double& mesonRadius,
			DoubleParameter& couplingHidden, DoubleParameter& coupling,
			double _massHiddenChannelA, double _massHiddenChannelB,
			int _subsys, int resSpin, int m, int n) ;

	AmpFlatteRes(const AmpFlatteRes&, const char*);
	AmpFlatteRes(const AmpFlatteRes&);

	virtual ~AmpFlatteRes();

	void setBarrierMass(double, double);

	virtual void initialise();
	virtual std::complex<double> evaluate() const;
	//  virtual double evaluate(double x[],int dim, void * param) const {return 0;};//used for MC integration
	void setDecayMasses(double, double, double, double);
	//  double getMaximum() const{return 1;};
	double integral() const {return 1;};
	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	inline virtual bool isSubSys(const unsigned int subSys)const{return (subSys==_subSys);};

protected:
	DoubleParameter _couplingHiddenChannel;
	DoubleParameter _coupling;
	AmpWigner _wignerD;

	double _massHiddenChannelA;//hidden channel: mass particle A
	double _massHiddenChannelB; //hidden channel: mass particle B

private:
	//ClassDef(AmpFlatteRes,1) // Relativistic Breit-Wigner resonance model

};

#endif
