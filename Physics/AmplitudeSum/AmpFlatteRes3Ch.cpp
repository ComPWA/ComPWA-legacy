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
		std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
		std::shared_ptr<DoubleParameter> mass, int subSys, Spin spin, Spin m, Spin n,
		std::shared_ptr<DoubleParameter> mesonRadius,
		std::shared_ptr<DoubleParameter> motherRadius,
		std::shared_ptr<DoubleParameter> g1,
		std::shared_ptr<DoubleParameter> g2, double g2_partA, double g2_partB,
		std::shared_ptr<DoubleParameter> g3, double g3_partA, double g3_partB,
		formFactorType type,
		int nCalls, normStyle nS) :
		AmpAbsDynamicalFunction(name, mag, phase, mass, subSys, spin, m, n,
				mesonRadius, motherRadius, type, nCalls, nS),
				_g1(g1),
				_g2(g2),_g2_partA(g2_partA), _g2_partB(g2_partB),
				_g3(g3),_g3_partA(g3_partA), _g3_partB(g3_partB)
{
	if(_g2_partA<0||_g2_partA>5||_g2_partB<0||_g2_partB>5)
		throw std::runtime_error("AmpFlatteRes3Ch::evaluateAmp | particle masses for second channel not set!");
	if(_g3_partA<0||_g3_partA>5||_g3_partB<0||_g3_partB>5)
		throw std::runtime_error("AmpFlatteRes3Ch::evaluateAmp | particle masses for third channel not set!");

	if( _g1->GetValue()!=tmp_g1 || _g2->GetValue()!=tmp_g2 || _g3->GetValue()!=tmp_g3 ) {
		SetModified();
		tmp_g1 = _g1->GetValue();
		tmp_g2 = _g2->GetValue();
		tmp_g3 = _g3->GetValue();
	}

	//setting default normalization
	GetNormalization();
}

AmpFlatteRes3Ch::~AmpFlatteRes3Ch()
{
}

//std::complex<double> AmpFlatteRes3Ch::evaluate(dataPoint& point) {
//	return AmpAbsDynamicalFunction::evaluate(point);
//}

std::complex<double> AmpFlatteRes3Ch::evaluateAmp(dataPoint& point) {
	if( _g1->GetValue()!=tmp_g1 || _g2->GetValue()!=tmp_g2 || _g3->GetValue()!=tmp_g3 ) {
		SetModified();
		tmp_g1 = _g1->GetValue();
		tmp_g2 = _g2->GetValue();
		tmp_g3 = _g3->GetValue();
	}

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double mSq=-999;
	switch(_subSys){
	case 3:
		mSq=(kin->getThirdVariableSq(point.getVal(0),point.getVal(1)));	break;
	case 4: mSq=(point.getVal(1)); break;
	case 5: mSq=(point.getVal(0)); break;
	}

	return dynamicalFunction(mSq,_mass->GetValue(),_ma,_mb,_g1->GetValue(),
			_g2_partA,_g2_partB,_g2->GetValue(),
			_g3_partA,_g3_partB,_g3->GetValue(),
			_spin,_mesonRadius->GetValue(), ffType);
}

std::complex<double> AmpFlatteRes3Ch::dynamicalFunction(double mSq, double mR,
		double massA1, double massA2, double gA,
		double massB1, double massB2, double couplingB,
		double massC1, double massC2, double couplingC,
		unsigned int J, double mesonRadius, formFactorType ffType ){
	std::complex<double> i(0,1);
	double sqrtS = sqrt(mSq);

	//channel A - signal channel
	std:complex<double> gammaA, qTermA, termA;
	double barrierA;
	//break-up momentum
	barrierA = Kinematics::FormFactor(sqrtS,massA1,massA2,J,mesonRadius, ffType)/
			Kinematics::FormFactor(mR,massA1,massA2,J,mesonRadius, ffType);
	//convert coupling to partial width of channel A
	gammaA = couplingToWidth(mSq,mR,gA,massA1,massA2,J,mesonRadius, ffType);
	//including the factor qTermA, as suggested by PDG, leads to an amplitude that doesn't converge.
	//		qTermA = Kinematics::qValue(sqrtS,massA1,massA2) / Kinematics::qValue(mR,massA1,massA2);
	qTermA = std::complex<double>(1,0);
	termA = gammaA*barrierA*barrierA*std::pow(qTermA,(double)2*J+1);

	//	std::cout<<qTermA<<" "<<Kinematics::qValue(sqrtS,massA1,massA2)
	//	<<" "<<Kinematics::qValue(mR,massA1,massA2)<<std::endl;

	//channel B - hidden channel
	std::complex<double> gammaB, qTermB, termB;
	double barrierB, gB;
	//break-up momentum
	barrierB = Kinematics::FormFactor(sqrtS,massB1,massB2,J,mesonRadius, ffType)/
			Kinematics::FormFactor(mR,massB1,massB2,J,mesonRadius, ffType);
	gB = couplingB;
	//convert coupling to partial width of channel B
	gammaB = couplingToWidth(mSq,mR,gB,massB1,massB2,J,mesonRadius, ffType);
	//		qTermB = Kinematics::qValue(sqrtS,massB1,massB2) / Kinematics::qValue(mR,massB1,massB2);
	qTermB = std::complex<double>(1,0);
	termB = gammaB*barrierB*barrierB*std::pow(qTermB,(double)2*J+1);

	//channel C - hidden channel
	std::complex<double> gammaC, qTermC, termC;
	double barrierC, gC;
	//break-up momentum
	barrierC = Kinematics::FormFactor(sqrtS,massC1,massC2,J,mesonRadius, ffType)/
			Kinematics::FormFactor(mR,massC1,massC2,J,mesonRadius, ffType);
	gC = couplingC;
	//convert coupling to partial width of channel C
	gammaC = couplingToWidth(mSq,mR,gC,massC1,massC2,J,mesonRadius, ffType);
	//		qTermC = Kinematics::qValue(sqrtS,massC1,massC2) / Kinematics::qValue(mR,massC1,massC2);
	qTermC = std::complex<double>(1,0);
	termC = gammaC*barrierC*barrierC*std::pow(qTermC,(double)2*J+1);

	//Coupling constant from production reaction. In case of a particle decay the production
	//coupling doesn't depend in energy since the CM energy is in the (RC) system fixed to the
	//mass of the decaying particle
	double g_production = 1;

	//-- old approach(BaBar)
	//std::complex<double> denom( mR*mR - mSq, (-1)*(rhoA*gA*gA + rhoB*gB*gB + rhoC*gC*gC) );
	//-- new approach - for spin 0 resonances in the imaginary part of the denominator the term qTerm
	//is added, compared to the old formula
	std::complex<double> denom = std::complex<double>( mR*mR - mSq,0)
			+ (-1.0)*i*sqrtS*(termA + termB + termC);

	std::complex<double> result = std::complex<double>(gA*g_production,0) / denom;

	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"mpFlatteRes3Ch::dynamicalFunction() | "<<barrierA<<" "<<mR<<" "<<mSq<<" "
				<<massA1<<" "<<massA2<<std::endl;
		return 0;
	}
	return result;
}
std::shared_ptr<FunctionTree> AmpFlatteRes3Ch::setupTree(
		allMasses& theMasses,allMasses& toyPhspSample,std::string suffix, ParameterList& params){
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double phspVol = kin->getPhspVolume();
	BOOST_LOG_TRIVIAL(info) << "AmpFlatteRes3Ch::setupBasicTree() | "<<_name;
	//------------Setup Tree---------------------
	std::shared_ptr<FunctionTree> newTree(new FunctionTree());
	//------------Setup Tree Pars---------------------
	std::shared_ptr<MultiDouble> m23sq( new MultiDouble("m23sq",theMasses.masses_sq.at( std::make_pair(2,3) )) );
	std::shared_ptr<MultiDouble> m13sq( new MultiDouble("m13sq",theMasses.masses_sq.at( std::make_pair(1,3) )) );
	std::shared_ptr<MultiDouble> m12sq( new MultiDouble("m12sq",theMasses.masses_sq.at( std::make_pair(1,2) )) );
	std::shared_ptr<MultiDouble> m23sq_phsp( new MultiDouble("m23sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(2,3) )) );
	std::shared_ptr<MultiDouble> m13sq_phsp( new MultiDouble("m13sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(1,3) )) );
	std::shared_ptr<MultiDouble> m12sq_phsp( new MultiDouble("m12sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(1,2) )) );

	//----Strategies needed
	std::shared_ptr<MultAll> mmultStrat(new MultAll(ParType::MCOMPLEX));
	std::shared_ptr<MultAll> mmultDStrat(new MultAll(ParType::MDOUBLE));
	std::shared_ptr<AddAll> maddStrat(new AddAll(ParType::MCOMPLEX));
	std::shared_ptr<AbsSquare> msqStrat(new AbsSquare(ParType::MDOUBLE));
	std::shared_ptr<LogOf> mlogStrat(new LogOf(ParType::MDOUBLE));
	std::shared_ptr<MultAll> multStrat(new MultAll(ParType::COMPLEX));
	std::shared_ptr<MultAll> multDStrat(new MultAll(ParType::DOUBLE));
	std::shared_ptr<AddAll> addStrat(new AddAll(ParType::DOUBLE));
	std::shared_ptr<AddAll> addComplexStrat(new AddAll(ParType::COMPLEX));
	std::shared_ptr<AbsSquare> sqStrat(new AbsSquare(ParType::DOUBLE));
	std::shared_ptr<LogOf> logStrat(new LogOf(ParType::DOUBLE));
	std::shared_ptr<Complexify> complStrat(new Complexify(ParType::COMPLEX));
	std::shared_ptr<Inverse> invStrat(new Inverse(ParType::DOUBLE));
	std::shared_ptr<SquareRoot> sqRootStrat(new SquareRoot(ParType::DOUBLE));

	//----Add Nodes
	std::shared_ptr<WignerDStrategy> angdStrat(	new WignerDStrategy(_name,ParType::MDOUBLE) );
	std::shared_ptr<WignerDPhspStrategy> angdPhspStrat(	new WignerDPhspStrategy(_name,ParType::MDOUBLE) );

	std::shared_ptr<Flatte3ChStrategy> flatteStrat(new Flatte3ChStrategy(_name,ParType::MCOMPLEX));
	std::shared_ptr<Flatte3ChPhspStrategy> flattePhspStrat(new Flatte3ChPhspStrategy(_name,ParType::MCOMPLEX));

	newTree->createHead("Reso_"+_name, mmultStrat, theMasses.nEvents); //Reso=BW*C_*AD*N_

	newTree->createNode("FlatteRes_"+_name, flatteStrat, "Reso_"+_name, theMasses.nEvents); //BW
	newTree->createNode("C_"+_name, complStrat, "Reso_"+_name); //c
	newTree->createLeaf("Intens_"+_name, params.GetDoubleParameter("mag_"+_name), "C_"+_name); //r
	newTree->createLeaf("Phase_"+_name, params.GetDoubleParameter("phase_"+_name), "C_"+_name); //phi
	newTree->createNode("AngD_"+_name, angdStrat, "Reso_"+_name, theMasses.nEvents); //AD

	//Flatte
	newTree->createLeaf("m23sq", m23sq, "FlatteRes_"+_name); //ma
	newTree->createLeaf("m13sq", m13sq, "FlatteRes_"+_name); //mb
	newTree->createLeaf("m12sq", m12sq, "FlatteRes_"+_name); //mc
	newTree->createLeaf("subSysFlag_"+_name, _subSys, "FlatteRes_"+_name); //subSysFlag
	newTree->createLeaf("spin_"+_name, _spin, "FlatteRes_"+_name); //spin
	newTree->createLeaf("d_"+_name, params.GetDoubleParameter("d_"+_name) , "FlatteRes_"+_name); //d
	newTree->createLeaf("ffType_"+_name, ffType , "FlatteRes_"+_name); //d
	if(_name.find("a_0(980)") != _name.npos){
		try {
			newTree->createLeaf("m0_a_0", params.GetDoubleParameter("m0_a_0"), "FlatteRes_"+_name);//use global parameter g1_a0 (asdfef)
		} catch (BadParameter& e) {
			newTree->createLeaf("m0_"+_name, params.GetDoubleParameter("m0_"+_name), "FlatteRes_"+_name);//use local parameter g1_a0
		}
		try {
			newTree->createLeaf("g1_a_0", params.GetDoubleParameter("g1_a_0"), "FlatteRes_"+_name);//use global parameter g1_a0 (asdfef)
		} catch (BadParameter& e) {
			newTree->createLeaf("g1_"+_name, params.GetDoubleParameter("g1_"+_name), "FlatteRes_"+_name);//use local parameter g1_a0
		}
		try {
			newTree->createLeaf("g2_a_0", params.GetDoubleParameter("g2_a_0"), "FlatteRes_"+_name);//use global parameter g1_a0 (asdfef)
		} catch (BadParameter& e) {
			newTree->createLeaf("g2_"+_name, params.GetDoubleParameter("g2_"+_name), "FlatteRes_"+_name);
		}
		try {
			newTree->createLeaf("g3_a_0", params.GetDoubleParameter("g1_a_0"), "FlatteRes_"+_name);//use global parameter g1_a0 (asdfef)
		} catch (BadParameter& e) {
			newTree->createLeaf("g3_"+_name, params.GetDoubleParameter("g1_"+_name), "FlatteRes_"+_name);
		}
	} else {
		newTree->createLeaf("m0_"+_name, params.GetDoubleParameter("m0_"+_name), "FlatteRes_"+_name); //m0
		newTree->createLeaf("g1_"+_name, params.GetDoubleParameter("g1_"+_name), "FlatteRes_"+_name);
		newTree->createLeaf("g2_"+_name, params.GetDoubleParameter("g2_"+_name), "FlatteRes_"+_name);
		try {
			newTree->createLeaf("g3_"+_name, params.GetDoubleParameter("g3_"+_name), "FlatteRes_"+_name);
		} catch (BadParameter& e) {
			newTree->createLeaf("g3_"+_name, params.GetDoubleParameter("g1_"+_name), "FlatteRes_"+_name);
		}
	}
	newTree->createLeaf("massB1_"+_name, _g2_partA, "FlatteRes_"+_name);
	newTree->createLeaf("massB2_"+_name, _g2_partB, "FlatteRes_"+_name);
	newTree->createLeaf("massC1_"+_name, _g3_partA, "FlatteRes_"+_name);
	newTree->createLeaf("massC2_"+_name, _g3_partB, "FlatteRes_"+_name);
	//Angular distribution
	newTree->createLeaf("m23sq", m23sq, "AngD_"+_name); //ma
	newTree->createLeaf("m13sq", m13sq, "AngD_"+_name); //mb
	newTree->createLeaf("m12sq", m12sq, "AngD_"+_name); //mc
	newTree->createLeaf("M", kin->M, "AngD_"+_name); //M
	newTree->createLeaf("m1", kin->m1, "AngD_"+_name); //m1
	newTree->createLeaf("m2", kin->m2, "AngD_"+_name); //m2
	newTree->createLeaf("m3", kin->m3, "AngD_"+_name); //m3
	newTree->createLeaf("subSysFlag_"+_name, _subSys, "AngD_"+_name); //subSysFlag
	newTree->createLeaf("spin_"+_name,_spin, "AngD_"+_name); //spin
	newTree->createLeaf("m_"+_name, 0, "AngD_"+_name); //OutSpin 1
	newTree->createLeaf("n_"+_name, 0, "AngD_"+_name); //OutSpin 2

	//Normalization
	if(_normStyle!=normStyle::none){
		newTree->createNode("N_"+_name, sqRootStrat, "Reso_"+_name); //N = sqrt(NSq)
		newTree->createNode("NSq_"+_name, multDStrat, "N_"+_name); //NSq = N_phspMC * 1/PhspVolume * 1/Sum(|A|^2)
		newTree->createLeaf("PhspSize_"+_name, toyPhspSample.nEvents, "NSq_"+_name); // N_phspMC
		newTree->createLeaf("PhspVolume_"+_name, 1/phspVol, "NSq_"+_name); // 1/PhspVolume
		newTree->createNode("InvSum_"+_name, invStrat, "NSq_"+_name); //1/Sum(|A|^2)
		newTree->createNode("Sum_"+_name, addStrat, "InvSum_"+_name); //Sum(|A|^2)
		newTree->createNode("AbsVal_"+_name, msqStrat, "Sum_"+_name); //|A_i|^2
		newTree->createNode("NormReso_"+_name, mmultStrat, "AbsVal_"+_name, toyPhspSample.nEvents); //BW

		//Flatte (Normalization)
		newTree->createNode("NormFlatte_"+_name, flattePhspStrat, "NormReso_"+_name, toyPhspSample.nEvents); //BW
		newTree->createLeaf("m23sq_phsp", m23sq_phsp, "NormFlatte_"+_name); //ma
		newTree->createLeaf("m13sq_phsp", m13sq_phsp, "NormFlatte_"+_name); //mb
		newTree->createLeaf("m12sq_phsp", m12sq_phsp, "NormFlatte_"+_name); //mc
		newTree->createLeaf("subSysFlag_"+_name, _subSys, "NormFlatte_"+_name); //subSysFlag
		newTree->createLeaf("spin_"+_name, _spin, "NormFlatte_"+_name); //spin
		newTree->createLeaf("d_"+_name,  params.GetDoubleParameter("d_"+_name), "NormFlatte_"+_name); //d
		newTree->createLeaf("ffType_"+_name, ffType , "NormFlatte_"+_name); //d
		if(_name.find("a_0(980)") != _name.npos){
			try {
				newTree->createLeaf("m0_a_0", params.GetDoubleParameter("m0_a_0"), "NormFlatte_"+_name);//use global parameter m0_a0 (asdfef)
			} catch (BadParameter& e) {
				newTree->createLeaf("m0_"+_name, params.GetDoubleParameter("m0_"+_name), "NormFlatte_"+_name);//use local parameter m0_a0
			}
			try {
				newTree->createLeaf("g1_a_0", params.GetDoubleParameter("g1_a_0"), "NormFlatte_"+_name);//use global parameter g1_a0 (asdfef)
			} catch (BadParameter& e) {
				newTree->createLeaf("g1_"+_name, params.GetDoubleParameter("g1_"+_name), "NormFlatte_"+_name);//use local parameter g1_a0
			}
			try {
				newTree->createLeaf("g2_a_0", params.GetDoubleParameter("g2_a_0"), "NormFlatte_"+_name);//use global parameter g1_a0 (asdfef)
			} catch (BadParameter& e) {
				newTree->createLeaf("g2_"+_name, params.GetDoubleParameter("g2_"+_name), "NormFlatte_"+_name);
			}
			try {
				newTree->createLeaf("g3_a_0", params.GetDoubleParameter("g1_a_0"), "NormFlatte_"+_name);//use global parameter g1_a0 (asdfef)
			} catch (BadParameter& e) {
				newTree->createLeaf("g3_"+_name, params.GetDoubleParameter("g1_"+_name), "NormFlatte_"+_name);
			}
		} else {
			newTree->createLeaf("m0_"+_name, params.GetDoubleParameter("m0_"+_name), "NormFlatte_"+_name); //m0
			newTree->createLeaf("g1_"+_name, params.GetDoubleParameter("g1_"+_name), "NormFlatte_"+_name);//use local parameter g1_a0
			newTree->createLeaf("g2_"+_name, params.GetDoubleParameter("g2_"+_name), "NormFlatte_"+_name);
			try {
				newTree->createLeaf("g3_"+_name, params.GetDoubleParameter("g3_"+_name), "NormFlatte_"+_name);
			} catch (BadParameter& e) {
				newTree->createLeaf("g3_"+_name, params.GetDoubleParameter("g1_"+_name), "NormFlatte_"+_name);
			}
		}

		newTree->createLeaf("massB1_"+_name, _g2_partA, "NormFlatte_"+_name);
		newTree->createLeaf("massB2_"+_name, _g2_partB, "NormFlatte_"+_name);
		newTree->createLeaf("massC1_"+_name, _g3_partA, "NormFlatte_"+_name);
		newTree->createLeaf("massC2_"+_name, _g3_partB, "NormFlatte_"+_name);

		//Angular distribution (Normalization)
		newTree->createNode("NormAngD_"+_name, angdPhspStrat, "NormReso_"+_name, toyPhspSample.nEvents); //AD
		newTree->createLeaf("m23sq_phsp", m23sq_phsp, "NormAngD_"+_name); //ma
		newTree->createLeaf("m13sq_phsp", m13sq_phsp, "NormAngD_"+_name); //mb
		newTree->createLeaf("m12sq_phsp", m12sq_phsp, "NormAngD_"+_name); //mc
		newTree->createLeaf("M", kin->M, "NormAngD_"+_name); //M
		newTree->createLeaf("m1", kin->m1, "NormAngD_"+_name); //m1
		newTree->createLeaf("m2", kin->m2, "NormAngD_"+_name); //m2
		newTree->createLeaf("m3", kin->m3, "NormAngD_"+_name); //m3
		newTree->createLeaf("subSysFlag_"+_name, _subSys, "NormAngD_"+_name); //subSysFlag
		newTree->createLeaf("spin_"+_name,_spin, "NormAngD_"+_name); //spin
		newTree->createLeaf("m_"+_name, 0, "NormAngD_"+_name); //OutSpin 1
		newTree->createLeaf("n_"+_name, 0, "NormAngD_"+_name); //OutSpin 2
	} else {
		newTree->createLeaf("N_"+_name, 1., "Reso_"+_name);
	}

	switch(_subSys){
	case 3:{ //reso in sys of particles 1&2
		//newTree->createLeaf("mym_"+_name, m12, "RelBW_"+_name); //m
		newTree->createLeaf("ma_"+_name, kin->m1, "FlatteRes_"+_name); //ma
		newTree->createLeaf("mb_"+_name, kin->m2, "FlatteRes_"+_name); //mb
		if(_normStyle!=normStyle::none){
			newTree->createLeaf("ma_"+_name, kin->m1, "NormFlatte_"+_name); //ma
			newTree->createLeaf("mb_"+_name, kin->m2, "NormFlatte_"+_name); //mb
		}
		break;
	}
	case 4:{ //reso in sys of particles 1&3
		//newTree->createLeaf("mym_"+_name, m13, "FlatteRes_"+_name); //m
		newTree->createLeaf("ma_"+_name, kin->m1, "FlatteRes_"+_name); //ma
		newTree->createLeaf("mb_"+_name, kin->m3, "FlatteRes_"+_name); //mb
		if(_normStyle!=normStyle::none){
			newTree->createLeaf("ma_"+_name, kin->m1, "NormFlatte_"+_name); //ma
			newTree->createLeaf("mb_"+_name, kin->m3, "NormFlatte_"+_name); //mb
		}
		break;
	}
	case 5:{ //reso in sys of particles 2&3
		//newTree->createLeaf("mym_"+_name, m23, "FlatteRes_"+_name); //m
		newTree->createLeaf("ma_"+_name, kin->m2, "FlatteRes_"+_name); //ma
		newTree->createLeaf("mb_"+_name, kin->m3, "FlatteRes_"+_name); //mb
		if(_normStyle!=normStyle::none){
			newTree->createLeaf("ma_"+_name, kin->m2, "NormFlatte_"+_name); //ma
			newTree->createLeaf("mb_"+_name, kin->m3, "NormFlatte_"+_name); //mb
		}
		break;
	}
	default:{
		BOOST_LOG_TRIVIAL(error)<<"AmpSumIntensity::setupBasicTree(): Subsys not found!!";
	}
	}
	return newTree;
}



bool Flatte3ChStrategy::execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
	if( checkType != out->type() ) {
		throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
		return false;
	}

	double m0, d,  ma, mb, g1, gB, massB1, massB2, gC, massC1, massC2;
	unsigned int spin, subSys;
	int ffType;
	//Get parameters from ParameterList -
	//enclosing in try...catch for the case that names of nodes have changed
	try{
		m0 = double(paras.GetParameterValue("m0_a_0"));//special case for peter's channel
	}catch(BadParameter& e){
		try{
			m0 = double(paras.GetParameterValue("m0_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter m0_a_0";
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter m0_"+name;
			throw;
		}
	}
	try{
		spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter ParOfNode_spin_"+name;
		throw;
	}
	try{
		d = double(paras.GetParameterValue("d_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter d_"+name;
		throw;
	}
	try{
		ffType = double(paras.GetParameterValue("ParOfNode_ffType_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter ffType_"+name;
		throw;
	}
	//		norm = double(paras.GetParameterValue("ParOfNode_norm_"+name));
	try{
		subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter ParOfNode_subSysFlag_"+name;
		throw;
	}

	try{
		ma = double(paras.GetParameterValue("ParOfNode_ma_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter ParOfNode_ma_"+name;
		throw;
	}
	try{
		mb = double(paras.GetParameterValue("ParOfNode_mb_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter ParOfNode_mb_"+name;
		throw;
	}
	try{
		g1 = double(paras.GetParameterValue("g1_a_0"));//special case for peter's channel
	}catch(BadParameter& e){
		try{
			g1 = double(paras.GetParameterValue("g1_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter g1_a_0";
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter g1_"+name;
			throw;
		}
	}
	try{
		massB1 = double(paras.GetParameterValue("ParOfNode_massB1_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter ParOfNode_massB1_"+name;
		throw;
	}
	try{
		massB2 = double(paras.GetParameterValue("ParOfNode_massB2_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter ParOfNode_massB2_"+name;
		throw;
	}
	try{
		gB = double(paras.GetParameterValue("g2_a_0"));//special case for peter's channel
	}catch(BadParameter& e){
		try{
			gB = double(paras.GetParameterValue("g2_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter g2_a_0";
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter g2_"+name;
			throw;
		}
	}
	try{
		massC1 = double(paras.GetParameterValue("ParOfNode_massC1_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter ParOfNode_massC1_"+name;
		throw;
	}
	try{
		massC2 = double(paras.GetParameterValue("ParOfNode_massC2_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter ParOfNode_massC2"+name;
		throw;
	}
	try{
		gC = double(paras.GetParameterValue("g3_"+name));
	}catch(BadParameter& e){
		try{
			gC = double(paras.GetParameterValue("g1_a_0"));
		}catch(BadParameter& e){
			try{
				gC = double(paras.GetParameterValue("g1_"+name));
			}catch(BadParameter& e){
				BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter g3_"+name;
				BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter g1_a_0";
				BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter g1_"+name;
				throw;
			}
		}
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
						spin,d, formFactorType(ffType));
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
	case 3:{ mSq  = (double(paras.GetParameterValue("m12sq"))); break; }
	case 4:{ mSq  = (double(paras.GetParameterValue("m13sq"))); break; }
	case 5:{ mSq  = (double(paras.GetParameterValue("m23sq"))); break; }
	}

	std::complex<double> result = AmpFlatteRes3Ch::dynamicalFunction(mSq,m0,ma,mb,g1,
			massB1,massB2,gB,
			massC1,massC2,gC,
			spin,d, formFactorType(ffType));
	out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
	return true;
}

bool Flatte3ChPhspStrategy::execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
	if( checkType != out->type() ) {
		throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
		return false;
	}

	double m0, d, ma, mb, g1, gB, massB1, massB2, gC, massC1, massC2;
	unsigned int spin, subSys;
	int ffType;

	//Get parameters from ParameterList -
	//enclosing in try...catch for the case that names of nodes have changed
	try{
		m0 = double(paras.GetParameterValue("m0_a_0"));//special case for peter's channel
	}catch(BadParameter& e){
		try{
			m0 = double(paras.GetParameterValue("m0_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter m0_a_0";
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter m0_"+name;
			throw;
		}
	}
	try{
		spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter ParOfNode_spin_"+name;
		throw;
	}
	try{
		d = double(paras.GetParameterValue("d_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter d_"+name;
		throw;
	}
	try{
		ffType = double(paras.GetParameterValue("ParOfNode_ffType_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter ffType_"+name;
		throw;
	}

	try{
		subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter ParOfNode_subSysFlag_"+name;
		throw;
	}
	try{
		ma = double(paras.GetParameterValue("ParOfNode_ma_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter ParOfNode_ma_"+name;
		throw;
	}
	try{
		mb = double(paras.GetParameterValue("ParOfNode_mb_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter ParOfNode_mb_"+name;
		throw;
	}
	try{
		g1 = double(paras.GetParameterValue("g1_a_0"));
	}catch(BadParameter& e){
		try{
			g1 = double(paras.GetParameterValue("g1_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter g1_a_0";
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChStrategy: can't find parameter g1_"+name;
			throw;
		}
	}
	try{
		massB1 = double(paras.GetParameterValue("ParOfNode_massB1_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter ParOfNode_massB1_"+name;
		throw;
	}
	try{
		massB2 = double(paras.GetParameterValue("ParOfNode_massB2_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter ParOfNode_massB2_"+name;
		throw;
	}
	try{
		gB = double(paras.GetParameterValue("g2_a_0"));//special case for peter's channel
	}catch(BadParameter& e){
		try{
			gB = double(paras.GetParameterValue("g2_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter g2_a_0";
			BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter g2_"+name;
			throw;
		}
	}
	try{
		massC1 = double(paras.GetParameterValue("ParOfNode_massC1_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter ParOfNode_massC1_"+name;
		throw;
	}
	try{
		massC2 = double(paras.GetParameterValue("ParOfNode_massC2_"+name));
	}catch(BadParameter& e){
		BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter ParOfNode_massC2_"+name;
		throw;
	}
	try{
		gC = double(paras.GetParameterValue("g3_"+name));
	}catch(BadParameter& e){
		try{
			gC = double(paras.GetParameterValue("g1_a_0"));
		}catch(BadParameter& e){
			try{
				gC = double(paras.GetParameterValue("g1_"+name));
			}catch(BadParameter& e){
				BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter g3_"+name;
				BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter g1_a_0";
				BOOST_LOG_TRIVIAL(error) <<"Flatte3ChPhspStrategy: can't find parameter g1_"+name;
				throw;
			}
		}
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
						spin,d, formFactorType(ffType));
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
			spin,d, formFactorType(ffType));
	out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
	return true;
}

Flatte3ChConf::Flatte3ChConf(const boost::property_tree::ptree &pt_) : FlatteConf(pt_){
	m_g3= pt_.get<double>("g3");
	m_g3_fix= pt_.get<bool>("g3_fix");
	m_g3_min= pt_.get<double>("g3_min");
	m_g3_max= pt_.get<double>("g3_max");
	m_g3_part1= pt_.get<std::string>("g3_part1");
	m_g3_part2= pt_.get<std::string>("g3_part2");
}
void Flatte3ChConf::put(boost::property_tree::ptree &pt_){
	FlatteConf::put(pt_);
	pt_.put("g3", m_g3);
	pt_.put("g3_fix", m_g3_fix);
	pt_.put("g3_min", m_g3_min);
	pt_.put("g3_max", m_g3_max);
	pt_.put("g3_part1", m_g3_part1);
	pt_.put("g3_part2", m_g3_part2);
}
void Flatte3ChConf::update(ParameterList par){
	if(!m_enable) return;
	FlatteConf::update(par);

	std::shared_ptr<DoubleParameter> p;
	if(m_name.find("a_0(980)") != m_name.npos) {
		try{// only update parameters if they are found in list
			p = par.GetDoubleParameter("g1_a_0");
		} catch (BadParameter b) {
			try{// only update parameters if they are found in list
				p = par.GetDoubleParameter("g1_"+m_name);
			} catch (BadParameter b) {
			}
		}
	} else {
		try{// only update parameters if they are found in list
			p = par.GetDoubleParameter("g3_"+m_name);
		} catch (BadParameter b) {
			try{// only update parameters if they are found in list
				p = par.GetDoubleParameter("g1_"+m_name);
			} catch (BadParameter b) {
			}
		}
	}
	m_g3 = p->GetValue();
	m_g3_fix = p->IsFixed();
	m_g3_min = p->GetMinValue();
	m_g3_max = p->GetMaxValue();

	return;
}
