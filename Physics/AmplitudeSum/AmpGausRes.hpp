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

#ifndef AMP_GAUS_RES
#define AMP_GAUS_RES

#include <vector>

#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"

using namespace std;

class AmpGausRes : public AmpAbsDynamicalFunction  {
public:

	AmpGausRes(const char *name,
			std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
			std::shared_ptr<DoubleParameter> mass, int subSys,
			std::shared_ptr<DoubleParameter> width,
			int nCalls=30000, normStyle nS=normStyle::one) ;

	AmpGausRes(const AmpGausRes&);

	~AmpGausRes();

	virtual std::complex<double> evaluate(dataPoint& point);
	virtual std::complex<double> evaluateAmp(dataPoint& point);
	virtual double evaluateWignerD(dataPoint& point) const { return 1; };

	inline virtual bool isSubSys(const unsigned int subSys)const{return (subSys==_subSys);};
	double getSpin(){return 0;};

protected:
	std::shared_ptr<DoubleParameter> _width;

};

#endif
