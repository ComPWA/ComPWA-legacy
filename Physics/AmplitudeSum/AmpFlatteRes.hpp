//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding correct couplings
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
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

using namespace std;

class AmpFlatteRes : public AmpAbsDynamicalFunction, public AmpKinematics {
public:
	AmpFlatteRes(const char *name,
			DoubleParameter& resMass, DoubleParameter& resWidth,
			double& mesonRadius,
			DoubleParameter& coupling, DoubleParameter& couplingHidden,
			double _massHiddenChannelA, double _massHiddenChannelB,
			int _subsys, int resSpin, int m, int n) ;

	AmpFlatteRes(const AmpFlatteRes&, const char*);
	AmpFlatteRes(const AmpFlatteRes&);

	virtual ~AmpFlatteRes();

	void setBarrierMass(double, double);

	//static function for dynamic part
	static std::complex<double> dynamicalFunction(double mSq, double mR, double ma, double mb, double coupling,
		double mHiddenA, double mHiddenB, double couplingHidden,unsigned int J);

	virtual void initialise();
	std::complex<double> evaluate(dataPoint& point) { return ( _norm*evaluateAmp(point)*evaluateWignerD(point) ); }
	virtual std::complex<double> evaluateAmp(dataPoint& point) ;
	virtual double evaluateWignerD(dataPoint& point) { return _wignerD.evaluate(point); };

	void setDecayMasses(double, double, double, double);
	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	inline virtual bool isSubSys(const unsigned int subSys)const{return (subSys==_subSys);};

protected:
	DoubleParameter _couplingHiddenChannel;
	DoubleParameter _coupling;
//	AmpWigner _wignerD;
	AmpWigner2 _wignerD;

	double _massHiddenChannelA;//hidden channel: mass particle A
	double _massHiddenChannelB; //hidden channel: mass particle B

private:
	//ClassDef(AmpFlatteRes,1) // Relativistic Breit-Wigner resonance model

};


class FlatteStrategy : public Strategy {
public:
	FlatteStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}

	virtual const std::string to_str() const { return ("flatte amplitude of "+name); }

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
		if( checkType != out->type() ) {
			throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
			return false;
		}

		double m0, d;
		unsigned int spin, subSys;
		try{
			m0 = double(paras.GetParameterValue("m0_"+name));
		}catch(BadParameter& e){
			m0 = double(paras.GetParameterValue("ParOfNode_m0_"+name));
		}
		spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
		d = double(paras.GetParameterValue("ParOfNode_d_"+name));
		//		norm = double(paras.GetParameterValue("ParOfNode_norm_"+name));
		subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));

		double ma = double(paras.GetParameterValue("ParOfNode_ma_"+name));
		double mb = double(paras.GetParameterValue("ParOfNode_mb_"+name));
		double mHiddenA = double(paras.GetParameterValue("ParOfNode_mHiddenA_"+name));
		double mHiddenB = double(paras.GetParameterValue("ParOfNode_mHiddenB_"+name));
		double coupling = double(paras.GetParameterValue("ParOfNode_coupling_"+name));
		double couplingHidden = double(paras.GetParameterValue("ParOfNode_couplingHidden_"+name));

		//MultiDim output, must have multidim Paras in input
		if(checkType == ParType::MCOMPLEX){
			if(paras.GetNMultiDouble()){
				unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();

				std::vector<std::complex<double> > results(nElements, std::complex<double>(0.));
				std::shared_ptr<MultiDouble> mp;//=paras.GetMultiDouble("mym_"+name);
				switch(subSys){
				case 3:{ mp  = (paras.GetMultiDouble("m12sq")); break; }
				case 4:{ mp  = (paras.GetMultiDouble("m13sq")); break; }
				case 5:{ mp  = (paras.GetMultiDouble("m23sq")); break; }
				}

				//calc BW for each point
				for(unsigned int ele=0; ele<nElements; ele++){
					double mSq = (mp->GetValue(ele));
					results[ele] = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,coupling,mHiddenA,mHiddenB,couplingHidden,spin);
					//					if(ele<10) std::cout<<"Strategy BWrel "<<results[ele]<<std::endl;
				}

				//std::vector<std::complex<double> > resultsTMP(nElements, std::complex<double>(1.));
				out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));
				return true;
			}else{ //end multidim para treatment
				throw(WrongParType("Input MultiDoubles missing in BW strat of "+name));
				return false;
			}
		}//end multicomplex output


		double mSq;// = sqrt(paras.GetParameterValue("mym_"+name));
		switch(subSys){
		case 3:{ mSq  = (double(paras.GetParameterValue("m12sq"))); break; }
		case 4:{ mSq  = (double(paras.GetParameterValue("m13sq"))); break; }
		case 5:{ mSq  = (double(paras.GetParameterValue("m23sq"))); break; }
		}

		std::complex<double> result = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,coupling,mHiddenA,mHiddenB,couplingHidden,spin);
		out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
		return true;
	}

protected:
	std::string name;
};

class FlattePhspStrategy : public Strategy {
public:
	FlattePhspStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}

	virtual const std::string to_str() const { return ("flatte amplitude of "+name); }

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
		if( checkType != out->type() ) {
			throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
			return false;
		}

		double m0, d;
		unsigned int spin, subSys;
		try{
			m0 = double(paras.GetParameterValue("m0_"+name));
		}catch(BadParameter& e){
			m0 = double(paras.GetParameterValue("ParOfNode_m0_"+name));
		}
		spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
		d = double(paras.GetParameterValue("ParOfNode_d_"+name));
		//		norm = double(paras.GetParameterValue("ParOfNode_norm_"+name));
		subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));

		//m  = double(paras.GetParameterValue("mym"));
		double ma = double(paras.GetParameterValue("ParOfNode_ma_"+name));
		double mb = double(paras.GetParameterValue("ParOfNode_mb_"+name));
		double mHiddenA = double(paras.GetParameterValue("ParOfNode_mHiddenA_"+name));
		double mHiddenB = double(paras.GetParameterValue("ParOfNode_mHiddenB_"+name));
		double coupling = double(paras.GetParameterValue("ParOfNode_coupling_"+name));
		double couplingHidden = double(paras.GetParameterValue("ParOfNode_couplingHidden_"+name));

		//MultiDim output, must have multidim Paras in input
		if(checkType == ParType::MCOMPLEX){
			if(paras.GetNMultiDouble()){
				unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();

				std::vector<std::complex<double> > results(nElements, std::complex<double>(0.));
				std::shared_ptr<MultiDouble> mp;//=paras.GetMultiDouble("mym_"+name);
				switch(subSys){
				case 3:{ mp  = (paras.GetMultiDouble("m12sq_phsp")); break; }
				case 4:{ mp  = (paras.GetMultiDouble("m13sq_phsp")); break; }
				case 5:{ mp  = (paras.GetMultiDouble("m23sq_phsp")); break; }
				}

				//calc BW for each point
				for(unsigned int ele=0; ele<nElements; ele++){
					double mSq = (mp->GetValue(ele));
					results[ele] = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,coupling,mHiddenA,mHiddenB,couplingHidden,spin);
					//					if(ele<10) std::cout<<"Strategy BWrel "<<results[ele]<<std::endl;
				}
				out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));
				return true;
			}else{ //end multidim para treatment
				throw(WrongParType("Input MultiDoubles missing in BW strat of "+name));
				return false;
			}
		}//end multicomplex output


		double mSq;// = sqrt(paras.GetParameterValue("mym_"+name));
		switch(subSys){
		case 3:{ mSq  = (double(paras.GetParameterValue("m12sq_phsp"))); break; }
		case 4:{ mSq  = (double(paras.GetParameterValue("m13sq_phsp"))); break; }
		case 5:{ mSq  = (double(paras.GetParameterValue("m23sq_phsp"))); break; }
		}

		std::complex<double> result = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,coupling,mHiddenA,mHiddenB,couplingHidden,spin);
		out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
		return true;
	}

protected:
	std::string name;
};
#endif
