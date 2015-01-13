//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - adding correct g1s
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
			std::shared_ptr<DoubleParameter> resMass,
			std::shared_ptr<DoubleParameter> mesonRadius,
			std::shared_ptr<DoubleParameter> motherRadius,
			std::shared_ptr<DoubleParameter> g1, std::shared_ptr<DoubleParameter> g2,
			double _g2_partA, double _g2_partB,
			int _subsys, int resSpin, int m, int n) ;

	AmpFlatteRes(const AmpFlatteRes&, const char*);
	AmpFlatteRes(const AmpFlatteRes&);

	virtual ~AmpFlatteRes();

	void setBarrierMass(double, double);

	//static function for dynamic part
	static std::complex<double> dynamicalFunction(double mSq, double mR, double ma, double mb, double g1,
			double mHiddenA, double mHiddenB, double g2,unsigned int J);

	virtual void initialise() { };
	std::complex<double> evaluate(dataPoint& point) { return ( _norm*evaluateAmp(point)*evaluateWignerD(point) ); }
	virtual std::complex<double> evaluateAmp(dataPoint& point) ;
	virtual double evaluateWignerD(dataPoint& point) {
		if(_spin==0) return 1.0;//save some computing time
		return _wignerD.evaluate(point);
	};

	void setDecayMasses(double, double, double, double);
	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	inline virtual bool isSubSys(const unsigned int subSys)const{return (subSys==_subSys);};
	unsigned int getNParams(){ return nParams;}

protected:
	AmpWigner2 _wignerD;
	unsigned int nParams;

	double _g2_partA;//hidden channel: mass particle A
	double _g2_partB; //hidden channel: mass particle B
	std::shared_ptr<DoubleParameter> _g2, _g1;
	bool foundMasses;
	unsigned int id23, id13;
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

		double m0, d, ma, mb, g1, g2, mHiddenA, mHiddenB;
		unsigned int spin, subSys;
		//Get parameters from ParameterList -
		//enclosing in try...catch for the case that names of nodes have changed
		try{
			m0 = double(paras.GetParameterValue("m0_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter m0_"+name;
			throw;
		}
		try{
			spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_spin_"+name;
			throw;
		}
		try{
			d = double(paras.GetParameterValue("d_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter d_"+name;
			throw;
		}
		//		norm = double(paras.GetParameterValue("ParOfNode_norm_"+name));
		try{
			subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_subSysFlag_"+name;
			throw;
		}

		try{
			ma = double(paras.GetParameterValue("ParOfNode_ma_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_ma_"+name;
			throw;
		}

		try{
			mb = double(paras.GetParameterValue("ParOfNode_mb_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_mb_"+name;
			throw;
		}
		try{
			mHiddenA = double(paras.GetParameterValue("ParOfNode_mHiddenA_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_mHiddenA_"+name;
			throw;
		}
		try{
			mHiddenB = double(paras.GetParameterValue("ParOfNode_mHiddenB_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_mHiddenB_"+name;
			throw;
		}
		try{
			g1 = double(paras.GetParameterValue("g1_"+name));
		}catch(BadParameter& e){
			try{
				g1 = double(paras.GetParameterValue("g1_a_0"));//special case for peter's channel
			}catch(BadParameter& e){
				BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter g1_"+name;
				BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter g1_a_0";
				throw;
			}
		}
		try{
			g2 = double(paras.GetParameterValue("g2_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter g2_"+name;
			throw;
		}

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
					results[ele] = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,g1,mHiddenA,mHiddenB,g2,spin);
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

		std::complex<double> result = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,g1,mHiddenA,mHiddenB,g2,spin);
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

		double m0, d, ma, mb, g1, g2, mHiddenA, mHiddenB;
		unsigned int spin, subSys;

		//Get parameters from ParameterList -
		//enclosing in try...catch for the case that names of nodes have changed
		try{
			m0 = double(paras.GetParameterValue("m0_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter m0_"+name;
			throw;
		}
		try{
			spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter ParOfNode_spin_"+name;
			throw;
		}
		try{
			d = double(paras.GetParameterValue("d_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter d_"+name;
			throw;
		}
		try{
			subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter ParOfNode_subSysFlag_"+name;
			throw;
		}
		try{
			ma = double(paras.GetParameterValue("ParOfNode_ma_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter ParOfNode_ma_"+name;
			throw;
		}
		try{
			mb = double(paras.GetParameterValue("ParOfNode_mb_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter ParOfNode_mb_"+name;
			throw;
		}
		try{
			mHiddenA = double(paras.GetParameterValue("ParOfNode_mHiddenA_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter ParOfNode_mHiddenA_"+name;
			throw;
		}
		try{
			mHiddenB = double(paras.GetParameterValue("ParOfNode_mHiddenB_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter ParOfNode_mHiddenB_"+name;
			throw;
		}
		try{
			g1 = double(paras.GetParameterValue("g1_"+name));
		}catch(BadParameter& e){
			try{
				g1 = double(paras.GetParameterValue("g1_a_0"));
			}catch(BadParameter& e){
				BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter g1_"+name;
				BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter g1_a_0";
				throw;
			}
		}
		try{
			g2 = double(paras.GetParameterValue("g2_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter g2_"+name;
			throw;
		}


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
					results[ele] = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,g1,mHiddenA,mHiddenB,g2,spin);
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

		std::complex<double> result = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,g1,mHiddenA,mHiddenB,g2,spin);
		out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
		return true;
	}

protected:
	std::string name;
};
#endif
