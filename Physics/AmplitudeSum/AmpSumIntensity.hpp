//****************************************************************************
// Wrapper to provide intensity of amplitude sum
//****************************************************************************

#ifndef _AMPSUMINTENSITY_HPP
#define _AMPSUMINTENSITY_HPP

#include <vector>
#include <memory>

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooCmdArg.h"
#include "RooMsgService.h"
#include "RooGlobalFunc.h"
#include "RooCFunction1Binding.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooComplex.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooRealProxy.h"

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
	//constants
	//	const Double_t M; // GeV/c² (J/psi+)
	//	const Double_t Br; // GeV/c² (width)
	//	const Double_t m1; // GeV/c² (gamma)
	//	const Double_t m2; // GeV/c² (pi)
	//	const Double_t m3; // GeV/c² (pi)
	//	//const Double_t c = 299792458.; // m/s
	//	const Double_t PI; // m/s

	//	const Double_t m23_sq_min;
	//	const Double_t m23_sq_max;
	//	const Double_t m13_sq_min;
	//	const Double_t m13_sq_max;
	//	const Double_t m12_sq_min;
	//	const Double_t m12_sq_max;
	//
	//	const Double_t m23_min;
	//	const Double_t m23_max;
	//	const Double_t m13_min;
	//	const Double_t m13_max;
	//	const Double_t m12_min;
	//	const Double_t m12_max;

	RooRealVar ma;
	RooRealVar mb;
	RooRealVar mc;

	Double_t maxVal;

	AmpSumOfAmplitudes totAmp;

	//Resonance Variables
	std::vector<std::string> namer;
	std::vector<std::shared_ptr<RooRealVar> > mr;
	std::vector<std::shared_ptr<RooRealVar> > qr;
	std::vector<std::shared_ptr<RooRealVar> > gr;
	std::vector<std::shared_ptr<RooRealVar> > rr;
	std::vector<std::shared_ptr<RooRealVar> > par1;
	std::vector<std::shared_ptr<RooRealVar> > par2;
	std::vector<std::shared_ptr<RooRealVar> > phir;
	std::vector<std::shared_ptr<RooRealVar> > aj;
	std::vector<std::shared_ptr<RooRealVar> > am;
	std::vector<std::shared_ptr<RooRealVar> > an;

	std::vector<std::shared_ptr<AmpAbsDynamicalFunction> > rbw;

	std::vector<std::shared_ptr<AmpWigner> > angd;
	unsigned int nAmps;

	//	double lambda(double x, double y, double z){
	//		return x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z;
	//	}
	//
	//	double m13_sq_max_constr(double &m23_sq){
	//		return m1*m1+m3*m3+0.5/m23_sq*((M*M-m23_sq-m1*m1)*(m23_sq-m2*m2+m3*m3)+sqrt(lambda(m23_sq,M*M,m1*m1))*sqrt(lambda(m23_sq,m2*m2,m3*m3)));
	//	}
	//
	//	double m13_sq_min_constr(double &m23_sq){
	//		return m1*m1+m3*m3+0.5/m23_sq*((M*M-m23_sq-m1*m1)*(m23_sq-m2*m2+m3*m3)-sqrt(lambda(m23_sq,M*M,m1*m1))*sqrt(lambda(m23_sq,m2*m2,m3*m3)));
	//	}
	//
	//	double m23_sq_max_constr(double &m13_sq){
	//		return m2*m2+m3*m3+0.5/m13_sq*((M*M-m13_sq-m2*m2)*(m13_sq-m1*m1+m3*m3)+sqrt(lambda(m13_sq,M*M,m2*m2))*sqrt(lambda(m13_sq,m1*m1,m3*m3)));
	//	}
	//
	//	double m23_sq_min_constr(double &m13_sq){
	//		return m2*m2+m3*m3+0.5/m13_sq*((M*M-m13_sq-m2*m2)*(m13_sq-m1*m1+m3*m3)-sqrt(lambda(m13_sq,M*M,m2*m2))*sqrt(lambda(m13_sq,m1*m1,m3*m3)));
	//	}
	//
	//	double mb_sq_max_constr(double &ma_sq){
	//		return m13_sq_max_constr(ma_sq);
	//	}
	//
	//	double mb_sq_min_constr(double &ma_sq){
	//		return m13_sq_min_constr(ma_sq);
	//	}
	//
	//	double ma_sq_max_constr(double &mb_sq){
	//		return m23_sq_max_constr(mb_sq);
	//	}
	//
	//	double ma_sq_min_constr(double &mb_sq){
	//		return m23_sq_min_constr(mb_sq);
	//	}


private:


};

#endif
