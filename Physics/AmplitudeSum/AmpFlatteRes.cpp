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

#include <cmath>
#include <math.h>
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
AmpFlatteRes::AmpFlatteRes(const char *name,
		std::shared_ptr<DoubleParameter> resMass,
		std::shared_ptr<DoubleParameter> mesonRadius, //  meson radius
		std::shared_ptr<DoubleParameter> motherRadius, //  mother radius
		std::shared_ptr<DoubleParameter> g1, std::shared_ptr<DoubleParameter> g2,
		double g2_partA, double g2_partB, int subSys, int resSpin, int m, int n, double resRadius) :
		AmpAbsDynamicalFunction(name),
		AmpKinematics(resMass, subSys, resSpin, m, n, mesonRadius, motherRadius),
		_g2(g2), _g1(g1), _g2_partA(g2_partA), _g2_partB(g2_partB), mesonRadius(resRadius),
		foundMasses(false),	nParams(5)
{
	if(_g2_partA<0||_g2_partA>5||_g2_partB<0||_g2_partB>5)
		throw std::runtime_error("AmpFlatteRes3Ch::evaluateAmp | particle masses for second channel not set!");
}

AmpFlatteRes::~AmpFlatteRes() 
{
}

std::complex<double> AmpFlatteRes::evaluateAmp(dataPoint& point) {
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	if(!foundMasses){
		id23 = point.getID("m23sq");
		id13 = point.getID("m13sq");
		foundMasses=true;
	}
	double mSq=-999;
	switch(_subSys){
	case 3:
		mSq=(kin->getThirdVariableSq(point.getVal(0),point.getVal(1)));	break;
	case 4: mSq=(point.getVal(1)); break;
	case 5: mSq=(point.getVal(0)); break;
	}

	return dynamicalFunction(mSq,_mR->GetValue(),_ma,_mb,_g1->GetValue(), _g2_partA,_g2_partB,_g2->GetValue(),_spin,mesonRadius);
}
std::complex<double> AmpFlatteRes::dynamicalFunction(double mSq, double mR,
		double massA1, double massA2, double gA,
		double massB1, double massB2, double couplingRatio,
		unsigned int J, double mesonRadius ){
	double sqrtS = sqrt(mSq);

	//channel A - signal channel
	//break-up momentum
	double rhoA = AmpKinematics::phspFactor(sqrtS, massA1, massA2);
	double barrierA = AmpKinematics::FormFactor(sqrtS,massA1,massA2,J,mesonRadius)/AmpKinematics::FormFactor(mR,massA1,massA2,J,mesonRadius);
	//std::complex<double> qTermA = std::pow((qValue(sqrtS,massA1,massA2) / qValue(mR,massA1,massA2)), (2.*J+ 1.));
	//convert coupling to partial width of channel A
	double gammaA = couplingToWidth(mSq,mR,gA,massA1,massA2,J,mesonRadius);
	//including the factor qTermA, as suggested by PDG, leads to an amplitude that doesn't converge.
	//std::complex<double> termA = gammaA*qTermA*barrierA*barrierA;
	double termA = gammaA*barrierA*barrierA;

	//channel B - hidden channel
	//break-up momentum
	double rhoB = AmpKinematics::phspFactor(sqrtS, massB1, massB2);
	double barrierB = AmpKinematics::FormFactor(sqrtS,massB1,massB2,J,1.5)/AmpKinematics::FormFactor(mR,massB1,massB2,J,1.5);
	//std::complex<double> qTermB = std::pow((qValue(sqrtS,massB1,massB2) / qValue(mR,massB1,massB2)), (2.*J+ 1.));
	//we expect that the ratio (gA/gB)^2
	//double gB = sqrt(gA*gA/couplingRatio);
	double gB = couplingRatio; 
	//convert coupling to partial width of channel B
	double gammaB = couplingToWidth(mSq,mR,gB,massB1,massB2,J,mesonRadius);
	//std::complex<double> termB = gammaB*qTermB*barrierB*barrierB;
	double termB = gammaB*barrierB*barrierB;

	//Coupling constant from production reaction. In case of a particle decay the production
	//coupling doesn't depend in energy since the CM energy is in the (RC) system fixed to the
	//mass of the decaying particle
	double g_production = 1;

	//-- old approach(BaBar)
	//std::complex<double> denom( mR*mR - mSq, (-1)*(rhoA*gA*gA + rhoB*gB*gB) );
	//std::complex<double> result = std::complex<double>(gA,0) / denom;
	//-- new approach - for spin 0 resonances in the imaginary part of the denominator the term qTerm
	//is added, compared to the old formula
	std::complex<double> denom( mR*mR - mSq, (-1)*sqrtS*(termA + termB) );
	std::complex<double> result = std::complex<double>(gA*g_production,0) / denom;

	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"nan in Flatte: "<<barrierA<<" "<<mR<<" "<<mSq<<" "<<massA1<<" "<<massA2<<std::endl;
		return 0;
	}
	return result;
}

bool FlatteStrategy::execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
	if( checkType != out->type() ) {
		throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
		return false;
	}

	double m0, d, ma, mb, g1, g2, mHiddenA, mHiddenB;
	unsigned int spin, mesonRadius, subSys;
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
		mesonRadius = (double)(paras.GetParameterValue("ParOfNode_mesonRadius_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_mesonRadius_"+name;
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
				results[ele] = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,g1,mHiddenA,mHiddenB,g2,spin,mesonRadius);
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

	std::complex<double> result = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,g1,mHiddenA,mHiddenB,g2,spin,mesonRadius);
	out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
	return true;
}

bool FlattePhspStrategy::execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
	if( checkType != out->type() ) {
		throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
		return false;
	}

	double m0, d, ma, mb, g1, g2, mHiddenA, mHiddenB;
	unsigned int spin, subSys, mesonRadius;

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
		mesonRadius = (double)(paras.GetParameterValue("ParOfNode_mesonRadius_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_mesonRadius_"+name;
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
				results[ele] = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,g1,mHiddenA,mHiddenB,g2,spin,mesonRadius);
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

	std::complex<double> result = AmpFlatteRes::dynamicalFunction(mSq,m0,ma,mb,g1,mHiddenA,mHiddenB,g2,spin,mesonRadius);
	out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
	return true;
}
