//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding flatte type resonance, removing root dependence
//-------------------------------------------------------------------------------
//****************************************************************************
// Wrapper to provide intensity of amplitude sum
//****************************************************************************

#ifndef _AMPSUMINTENSITY_HPP
#define _AMPSUMINTENSITY_HPP

#include <vector>
#include <memory>
#include <map>
#include <string>

#include "Physics/Amplitude.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
//#include "Estimator/AmpFcn.cpp"
//#include "Optimizer/Minuit2/MinuitIF.hpp"
//#include "Estimator/MinLogLH/MinLogLH.hpp"

#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Core/Efficiency.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Generator.hpp"

class AmpSumIntensity : public Amplitude {

public:
	enum normStyle {
		none, /*!< no normaliztion between Amplitudes. */
		one, /*!< all amplitudes are normalized to one. The normalization factor is \f$ 1/\sqrt(\int |A|^2)\f$ */
		entries /*!<all amplitudes are normalized to the number of entries in dalitz plot. The normalization factor is \f$ 1/\sqrt(entries/area * \int |A|^2)\f$*/
	};
	//! Default Constructor (0x0)
	AmpSumIntensity(const double inM, const double inBr, const double in1,const double in2, const double in3,
			std::string nameM, std::string name1,std::string name2,std::string name3,
			 AmplitudeSetup ini, unsigned int entries=9999,
			normStyle ns=none, double dpArea=-999);
	AmpSumIntensity(AmplitudeSetup ini, normStyle ns, std::shared_ptr<Efficiency> eff=std::shared_ptr<Efficiency>(new UnitEfficiency()),
			 unsigned int entries=9999, double dpArea=-999);
	AmpSumIntensity(AmplitudeSetup ini, std::shared_ptr<Efficiency> eff=std::shared_ptr<Efficiency>(new UnitEfficiency()),
			unsigned int entries=9999, double dpArea=-999);
	AmpSumIntensity(const AmpSumIntensity& other);

	double evaluate(double x[], size_t dim);
	//! normalization integral for parameters \par
	virtual const double integral(ParameterList& par);
	//! normalization integral
	virtual const double integral();
	//! maximum value of amplitude with parameters \par
	virtual double getMaxVal( std::shared_ptr<Generator> gen);
	virtual double getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen);
	virtual void calcMaxVal(ParameterList& par ,std::shared_ptr<Generator> gen);
	virtual void calcMaxVal( std::shared_ptr<Generator> gen);

    virtual std::complex<double> getFirstAmp(dataPoint& point, ParameterList& par);
	virtual std::complex<double> getFirstReso(dataPoint& point, ParameterList& par);
    virtual std::complex<double> getFirstBW(dataPoint& point, ParameterList& par);


	virtual void setNevents(unsigned int n) { _entries=n; };
	virtual unsigned int getNevents() { return _entries; };

	//! setting new parameterList
	virtual void setParameterList(ParameterList& par);
	//! fill ParameterList with copied shared_ptr
	void copyParameterList(ParameterList& par);
	//! evaluate total amplitude using parameters \par at phsp point \point
	virtual const ParameterList& intensity(dataPoint& point, ParameterList& par);
	//! evaluate total amplitude using current set of parametersat phsp point \point
	virtual const ParameterList& intensity(dataPoint& point);
	virtual const ParameterList& intensity(std::vector<double> point, ParameterList& par);

//	virtual std::shared_ptr<FunctionTree> functionTree(allMasses& theMasses);
	virtual std::shared_ptr<FunctionTree> functionTree(allMasses& theMasses, allMasses& thePHSPMasses);
    virtual std::shared_ptr<FunctionTree> phspTree(allMasses& thePHSPMasses);
    virtual void resetTree(){
        myTree = std::shared_ptr<FunctionTree>();
        myPhspTree = std::shared_ptr<FunctionTree>();
    }

	virtual const bool fillStartParVec(ParameterList& outPar);

	virtual void printAmps();
	virtual double getFraction(std::string name) { return totAmp.getFraction(name)/integral(); };
	virtual double getFraction(unsigned int id) { return totAmp.getFraction(id)/integral(); };
	virtual void printFractions();
	virtual double getIntValue(std::string var1, double min1, double max1, std::string var2, double min2=0, double max2=0);

	double normReso(std::shared_ptr<AmpAbsDynamicalFunction> amp);
	virtual ~AmpSumIntensity(){};
	virtual AmpSumIntensity* Clone(){
		return (new AmpSumIntensity(*this));
	}

protected:
	void init();
//	void setupTree(allMasses& theMasses, bool isPhspTree);
	void setupTree(allMasses& thePHSPMasses){
		allMasses dummyMass;
		setupTree(thePHSPMasses, dummyMass);
	}
	void setupTree(allMasses& theMasses, allMasses& thePHSPMasses);

	std::shared_ptr<Efficiency> eff_;
	bool _calcMaxFcnVal;
	bool _calcNorm;
	double _maxFcnVal;
	AmpSumOfAmplitudes totAmp;
	AmplitudeSetup ampSetup;
	std::shared_ptr<FunctionTree> myTree;
	std::shared_ptr<FunctionTree> myPhspTree;
//	std::shared_ptr<ParameterList> treePar;

	double maxVal;

	normStyle _normStyle;
	unsigned int _entries;
	double _dpArea;
	unsigned int nAmps;

	//Resonance Variables
	std::vector<std::string> namer;
	std::map<std::string,std::shared_ptr<DoubleParameter> > mr;
	std::map<std::string,std::shared_ptr<DoubleParameter> > gr;
	std::map<std::string,std::shared_ptr<DoubleParameter> > rr;
	std::map<std::string,std::shared_ptr<DoubleParameter> > phir;

//	std::vector<std::shared_ptr<DoubleParameter> > qr;

//	std::vector<std::shared_ptr<IntegerParameter> > aj;
//	std::vector<std::shared_ptr<IntegerParameter> > am;
//	std::vector<std::shared_ptr<IntegerParameter> > an;

//	std::vector<std::shared_ptr<DoubleParameter> > par1;
//	std::vector<std::shared_ptr<DoubleParameter> > par2;

//	std::vector<std::shared_ptr<AmpWigner> > angd;

private:


};

#endif
