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
// Wrapper to provide intensity of amplitude sum
//****************************************************************************

#ifndef _AMPSUMINTENSITY_HPP
#define _AMPSUMINTENSITY_HPP

#include <vector>
#include <memory>

#include "Physics/Amplitude.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
//#include "Estimator/AmpFcn.cpp"
//#include "Optimizer/Minuit2/MinuitIF.hpp"
//#include "Estimator/MinLogLH/MinLogLH.hpp"

#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"
#include "Physics/DPKinematics/DPKinematics.hpp"
#include "Physics/DPKinematics/DataPoint.hpp"

class AmpSumIntensity : public Amplitude {

public:
	/// Default Constructor (0x0)
	AmpSumIntensity(const double inM, const double inBr, const double in1
			,const double in2, const double in3, AmplitudeSetup ini);
	AmpSumIntensity(DPKinematics kin, AmplitudeSetup ini);

	virtual double getMaxVal() { return 1; };

	virtual const double integral(ParameterList& par);

	virtual const ParameterList intensity(std::vector<double>& x, ParameterList& par);
	virtual const ParameterList intensity(ParameterList& par);

	virtual const bool fillStartParVec(ParameterList& outPar);

	virtual void printAmps();

	virtual ~AmpSumIntensity(){};

protected:
	void init();
	const DPKinematics _kin;
	AmpSumOfAmplitudes totAmp;
	AmplitudeSetup ampSetup;

	double maxVal;


	//Resonance Variables
	std::vector<std::string> namer;
	std::vector<std::shared_ptr<DoubleParameter> > mr;
	std::vector<std::shared_ptr<DoubleParameter> > gr;
	std::vector<std::shared_ptr<DoubleParameter> > rr;
	std::vector<std::shared_ptr<DoubleParameter> > phir;

	std::vector<std::shared_ptr<DoubleParameter> > qr;

	std::vector<std::shared_ptr<IntegerParameter> > aj;
	std::vector<std::shared_ptr<IntegerParameter> > am;
	std::vector<std::shared_ptr<IntegerParameter> > an;

	std::vector<std::shared_ptr<DoubleParameter> > par1;
	std::vector<std::shared_ptr<DoubleParameter> > par2;

	std::vector<std::shared_ptr<AmpWigner> > angd;
	unsigned int nAmps;

private:


};

#endif
