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
#include "Physics/AmplitudeSum/AmpFlatteRes3Ch.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"

AmpFlatteRes3Ch::AmpFlatteRes3Ch(const char *name,
		std::shared_ptr<DoubleParameter> resMass,
		std::shared_ptr<DoubleParameter> mesonRadius, //  meson radius
		std::shared_ptr<DoubleParameter> motherRadius, //  mother radius
		std::shared_ptr<DoubleParameter> g1,
		std::shared_ptr<DoubleParameter> g2, double g2_partA, double g2_partB,
		std::shared_ptr<DoubleParameter> g3, double g3_partA, double g3_partB,
		int subSys, int resSpin, int m, int n, double resRadius) :
		AmpAbsDynamicalFunction(name),
		AmpKinematics(resMass, subSys, resSpin, m, n, mesonRadius, motherRadius), _g1(g1),
		_g2(g2),_g2_partA(g2_partA), _g2_partB(g2_partB),
		_g3(g3),_g3_partA(g3_partA), _g3_partB(g3_partB),
		mesonRadius(resRadius),	foundMasses(false),	nParams(5)
{
	if(_g2_partA<0||_g2_partA>5||_g2_partB<0||_g2_partB>5)
		throw std::runtime_error("AmpFlatteRes3Ch::evaluateAmp | particle masses for second channel not set!");
	if(_g3_partA<0||_g3_partA>5||_g3_partB<0||_g3_partB>5)
		throw std::runtime_error("AmpFlatteRes3Ch::evaluateAmp | particle masses for third channel not set!");
}

AmpFlatteRes3Ch::~AmpFlatteRes3Ch()
{
}

std::complex<double> AmpFlatteRes3Ch::evaluateAmp(const dataPoint& point) {

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

	return dynamicalFunction(mSq,_mR->GetValue(),_ma,_mb,_g1->GetValue(),
			_g2_partA,_g2_partB,_g2->GetValue(),
			_g3_partA,_g3_partB,_g3->GetValue(),
			_spin,mesonRadius);
}
std::complex<double> AmpFlatteRes3Ch::dynamicalFunction(double mSq, double mR,
		double massA1, double massA2, double gA,
		double massB1, double massB2, double couplingRatioB,
		double massC1, double massC2, double couplingRatioC,
		unsigned int J, double mesonRadius ){
	std::complex<double> i(0,1);
	double sqrtS = sqrt(mSq);

	//channel A - signal channel
	//break-up momentum
//	std::complex<double> pA = AmpKinematics::phspFactor(sqrtS, massA1, massA2);
	double barrierA = AmpKinematics::FormFactor(sqrtS,massA1,massA2,J,mesonRadius)/AmpKinematics::FormFactor(mR,massA1,massA2,J,mesonRadius);
	//std::complex<double> qTermA = std::pow((qValue(sqrtS,massA1,massA2) / qValue(mR,massA1,massA2)), (2.*J+ 1.));
//	std::complex<double> qTermA = std::pow((phspFactor(sqrtS,massA1,massA2) / phspFactor(mR,massA1,massA2))*mR/sqrtS, (2*J+ 1));
	//convert coupling to partial width of channel A
	std::complex<double> gammaA = couplingToWidth(mSq,mR,gA,massA1,massA2,J,mesonRadius);
	//including the factor qTermA, as suggested by PDG, leads to an amplitude that doesn't converge.
	std::complex<double> termA = gammaA*barrierA*barrierA;

	//channel B - hidden channel
	//break-up momentum
//	std::complex<double> pB = AmpKinematics::phspFactor(sqrtS, massB1, massB2);
	double barrierB = AmpKinematics::FormFactor(sqrtS,massB1,massB2,J,1.5)/AmpKinematics::FormFactor(mR,massB1,massB2,J,1.5);
	//std::complex<double> qTermB = std::pow((qValue(sqrtS,massB1,massB2) / qValue(mR,massB1,massB2)), (2.*J+ 1.));
//	std::complex<double> qTermB = std::pow((phspFactor(sqrtS,massB1,massB2) / phspFactor(mR,massB1,massB2))*mR/sqrtS, (2*J+ 1));
	double gB = couplingRatioB;
	//	std::cout<<gA<< " "<<gB<<" "<<couplingRatio<<std::endl;
	//convert coupling to partial width of channel B
	std::complex<double> gammaB = couplingToWidth(mSq,mR,gB,massB1,massB2,J,mesonRadius);
	std::complex<double> termB = gammaB*barrierB*barrierB;

	//channel C - hidden channel
	//break-up momentum
//	std::complex<double> pC = AmpKinematics::phspFactor(sqrtS, massC1, massC2);
	double barrierC = AmpKinematics::FormFactor(sqrtS,massC1,massC2,J,1.5)/AmpKinematics::FormFactor(mR,massC1,massC2,J,1.5);
	//std::complex<double> qTermC = std::pow((qValue(sqrtS,massC1,massC2) / qValue(mR,massC1,massC2)), (2.*J+ 1.));
//	std::complex<double> qTermC = std::pow((phspFactor(sqrtS,massC1,massC2) / phspFactor(mR,massC1,massC2))*mR/sqrtS, (2*J+ 1));
	double gC = couplingRatioC;
	//convert coupling to partial width of channel C
	std::complex<double> gammaC = couplingToWidth(mSq,mR,gC,massC1,massC2,J,mesonRadius);
	std::complex<double> termC = gammaC*barrierC*barrierC;

	//Coupling constant from production reaction. In case of a particle decay the production
	//coupling doesn't depend in energy since the CM energy is in the (RC) system fixed to the
	//mass of the decaying particle
	double g_production = 1;

	//-- old approach(BaBar)
	//std::complex<double> denom( mR*mR - mSq, (-1)*(rhoA*gA*gA + rhoB*gB*gB + rhoC*gC*gC) );
	//-- new approach - for spin 0 resonances in the imaginary part of the denominator the term qTerm
	//is added, compared to the old formula
	std::complex<double> denom = std::complex<double>( mR*mR - mSq,0) + (-1.0)*i*sqrtS*(termA + termB+ termC);

	std::complex<double> result = std::complex<double>(gA*g_production,0) / denom;

	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"mpFlatteRes3Ch::dynamicalFunction() | "<<barrierA<<" "<<mR<<" "<<mSq<<" "
				<<massA1<<" "<<massA2<<std::endl;
		return 0;
	}
	return result;
}

bool Flatte3ChStrategy::execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
	if( checkType != out->type() ) {
		throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
		return false;
	}

	double m0, d, ma, mb, g1, gB, massB1, massB2, gC, massC1, massC2;
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
		massB1 = double(paras.GetParameterValue("ParOfNode_massB1_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_massB1_"+name;
		throw;
	}
	try{
		massB2 = double(paras.GetParameterValue("ParOfNode_massB2_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_massB2_"+name;
		throw;
	}
	try{
		gB = double(paras.GetParameterValue("g2_"+name));
	}catch(BadParameter& e){
		try{
			gB = double(paras.GetParameterValue("g1_a_0"));//special case for peter's channel
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter g1_"+name;
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter g1_a_0";
			throw;
		}
	}
	try{
		massC1 = double(paras.GetParameterValue("ParOfNode_massC1_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_massC1_"+name;
		throw;
	}
	try{
		massC2 = double(paras.GetParameterValue("ParOfNode_massC2_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter ParOfNode_massC2"+name;
		throw;
	}
	try{
		gC = double(paras.GetParameterValue("g3_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter g3_"+name;
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
				results[ele] = AmpFlatteRes3Ch::dynamicalFunction(mSq,m0,ma,mb,g1,
						massB1,massB2,gB,
						massC1,massC2,gC,
						spin,mesonRadius);
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

	std::complex<double> result = AmpFlatteRes3Ch::dynamicalFunction(mSq,m0,ma,mb,g1,
			massB1,massB2,gB,
			massC1,massC2,gC,
			spin,mesonRadius);
	out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
	return true;
}

bool Flatte3ChPhspStrategy::execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
	if( checkType != out->type() ) {
		throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
		return false;
	}

	double m0, d, ma, mb, g1, gB, massB1, massB2, gC, massC1, massC2;
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
		massB1 = double(paras.GetParameterValue("ParOfNode_massB1_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter ParOfNode_massB1_"+name;
		throw;
	}
	try{
		massB2 = double(paras.GetParameterValue("ParOfNode_massB2_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter ParOfNode_massB2_"+name;
		throw;
	}
	try{
		gB = double(paras.GetParameterValue("g2_"+name));
	}catch(BadParameter& e){
		try{
			gB = double(paras.GetParameterValue("g1_a_0"));//special case for peter's channel
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter g1_"+name;
			BOOST_LOG_TRIVIAL(error) <<"FlatteStrategy: can't find parameter g1_a_0";
			throw;
		}
	}
	try{
		massC1 = double(paras.GetParameterValue("ParOfNode_massC1_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter ParOfNode_massC1_"+name;
		throw;
	}
	try{
		massC2 = double(paras.GetParameterValue("ParOfNode_massC2_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter ParOfNode_massC2_"+name;
		throw;
	}
	try{
		gC = double(paras.GetParameterValue("g3_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"FlattePhspStrategy: can't find parameter g3_"+name;
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
				results[ele] = AmpFlatteRes3Ch::dynamicalFunction(mSq,m0,ma,mb,g1,
						massB1,massB2,gB,
						massC1,massC2,gC,
						spin,mesonRadius);
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

	std::complex<double> result = AmpFlatteRes3Ch::dynamicalFunction(mSq,m0,ma,mb,g1,
			massB1,massB2,gB,
			massC1,massC2,gC,
			spin,mesonRadius);
	out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
	return true;
}

Flatte3ChConf::Flatte3ChConf(const boost::property_tree::ptree &pt_) : FlatteConf(pt_){
	m_g3= pt_.get<double>("g3");
	m_g3_part1= pt_.get<std::string>("g3_part1");
	m_g3_part2= pt_.get<std::string>("g3_part2");
}
void Flatte3ChConf::put(boost::property_tree::ptree &pt_){
	FlatteConf::put(pt_);
	pt_.put("g3", m_g3);
	pt_.put("g3_part1", m_g3_part1);
	pt_.put("g3_part2", m_g3_part2);
}
void Flatte3ChConf::update(ParameterList par){
	FlatteConf::update(par);
	try{// only update parameters if they are found in list
	  if(m_name.find("a_0(980)") == 0)
		m_g2= par.GetDoubleParameter("g1_a_0")->GetValue();
	  else
		m_g2= par.GetDoubleParameter("g2_"+m_name)->GetValue();
	} catch (BadParameter b) { //do nothing if parameter is not found

	}
	try{// only update parameters if they are found in list
		m_g3= par.GetDoubleParameter("g3_"+m_name)->GetValue();
	} catch (BadParameter b) {
//		BOOST_LOG_TRIVIAL(error) <<"FlatteConf::update() | coupling g3 not found in parameter list!";
//		throw;
	}
}
