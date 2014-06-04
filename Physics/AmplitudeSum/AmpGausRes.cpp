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

std::complex<double> AmpGausRes::evaluateAmp(dataPoint& point) const {


	double m0 = _mR.GetValue();
	double width = _resWidth.GetValue();
	//  double m  = Double_t(_x);
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
//	double m23sq = point.getVal("m23sq");
//	double m13sq = point.getVal("m13sq");
	double m23sq = point.getVal(0);
	double m13sq = point.getVal(1);
	double m12sq = kin->getThirdVariableSq(m23sq,m13sq);
	double m = -999;
	switch(_subSys){
	case 3: m=sqrt(m12sq); break;
	case 4: m=sqrt(m13sq); break;
	case 5: m=sqrt(m23sq); break;
	}

	std::complex<double> gaus (_norm * exp(-1*(m-m0)*(m-m0)/width/width/2.),0);

	return gaus;
}
std::complex<double> AmpGausRes::evaluate(dataPoint& point) const {
	return evaluateAmp(point)*evaluateWignerD(point);
}
