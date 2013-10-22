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
#include "Estimator/AmpFcn.cpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"

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

	double getMaxVal();

	virtual const double integral(ParameterList& par);

	virtual const double intensity(std::vector<double>& x, ParameterList& par);
	virtual const double intensity(ParameterList& par);

	virtual const bool fillStartParVec(ParameterList& outPar);

	virtual void printAmps();

	virtual ~AmpSumIntensity(){};

protected:
	void init();
	AmplitudeSetup ampSetup;
	const DPKinematics _kin;

	double maxVal;

	AmpSumOfAmplitudes totAmp;

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
