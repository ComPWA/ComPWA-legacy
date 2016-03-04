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

AmpFlatteRes::AmpFlatteRes( normStyle nS, int calls ) :
		AmpAbsDynamicalFunction( nS, calls )
{

}

AmpFlatteRes::AmpFlatteRes(const char *name,
		std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
		std::shared_ptr<DoubleParameter> mass, int part1, int part2,
		Spin spin, Spin m, Spin n,
		std::shared_ptr<DoubleParameter> mesonRadius,
		std::shared_ptr<DoubleParameter> motherRadius,
		std::shared_ptr<DoubleParameter> g1,
		std::shared_ptr<DoubleParameter> g2, std::string g2_idA, std::string g2_idB,
		std::shared_ptr<DoubleParameter> g3, std::string g3_idA, std::string g3_idB,
		formFactorType type,
		int nCalls, normStyle nS) :
		AmpAbsDynamicalFunction(name, mag, phase, mass, part1, part2, spin, m, n,
				mesonRadius, motherRadius, type, nCalls, nS),
				_g1(g1),
				_g2(g2),_g2_idA(g2_idA), _g2_idB(g2_idB),
				_g3(g3),_g3_idA(g3_idA), _g3_idB(g3_idB)
{
	_g2_massA = PhysConst::instance()->getMass(_g2_idA);
	_g2_massB = PhysConst::instance()->getMass(_g2_idB);
	_g3_massA = PhysConst::instance()->getMass(_g3_idA);
	_g3_massB = PhysConst::instance()->getMass(_g3_idB);

	if(_g2_massA<0||_g2_massA>5||_g2_massB<0||_g2_massB>5)
		throw std::runtime_error("AmpFlatteRes::evaluateAmp | "
				"particle masses for second channel not set!");
	if(_g3_massA<0||_g3_massA>5||_g3_massB<0||_g3_massB>5)
		throw std::runtime_error("AmpFlatteRes::evaluateAmp | "
				"particle masses for third channel not set!");

	tmp_g1 = _g1->GetValue();
	tmp_g2 = _g2->GetValue();
	tmp_g3 = _g3->GetValue();

	initialize();
}

AmpFlatteRes::~AmpFlatteRes()
{
}

void AmpFlatteRes::Configure(boost::property_tree::ptree::value_type const& v,
		ParameterList& list)
{
	if( v.first != "Flatte" ) throw BadConfig("");

	boost::property_tree::ptree pt = v.second;
	AmpAbsDynamicalFunction::Configure(v,list);

	//Coupling to channel 1 (mandatory)
	auto tmp_g1_fix = pt.get<bool>("g1_fix",1);
	auto tmp_g1_min = pt.get<double>("g1_min",0.);
	auto tmp_g1_max = pt.get<double>("g1_max",10.);
	auto tmp_g1_name = pt.get_optional<std::string>("g1_name");
	if(!tmp_g1_name){
		auto tmp_g1= pt.get_optional<double>("g1");
		if(!tmp_g1)
			throw BadParameter("AmpFlatteRes::Configure() | "
					"g1 for "+_name+" not specified!");
		_g1 = std::shared_ptr<DoubleParameter>(
				new DoubleParameter(
						"g1_"+_name,tmp_g1.get(), tmp_g1_min, tmp_g1_max
				)
		);
		_g1->FixParameter(tmp_g1_fix);
		if(_enable) list.AddParameter(_g1);
		_g1_writeByName = 0;
	} else {
		try{
			_g1 = list.GetDoubleParameter(tmp_g1_name.get());
			_g1_writeByName = 1;
		} catch (BadParameter& ex){
			if(!_enable){
				_g1 = std::shared_ptr<DoubleParameter>(
						new DoubleParameter(
								tmp_g1_name.get(),0.0)
				);
				_g1_writeByName = 1;
			} else {
				BOOST_LOG_TRIVIAL(error) <<"AmpFlatteRes::Configure() | "
						"Requesting parameter "<<tmp_g1_name.get()<<" but"
						" was not found in parameter list. "
						"Quit since parameter is mandatory!";
				throw;
			}
		}
	}
	//Coupling to channel 2 (mandatory)
	auto tmp_g2_fix = pt.get<bool>("g2_fix",1);
	auto tmp_g2_min = pt.get<double>("g2_min",0.);
	auto tmp_g2_max = pt.get<double>("g2_max",10.);
	auto tmp_g2_name = pt.get_optional<std::string>("g2_name");
	if(!tmp_g2_name){
		auto tmp_g2= pt.get_optional<double>("g2");
		if(!tmp_g2)
			throw BadParameter("AmpFlatteRes::Configure() | "
					"g2 for "+_name+" not specified!");
		_g2 = std::shared_ptr<DoubleParameter>(
				new DoubleParameter(
						"g2_"+_name,tmp_g2.get(), tmp_g2_min, tmp_g2_max
				)
		);
		_g2->FixParameter(tmp_g2_fix);
		if(_enable) list.AddParameter(_g2);
		_g2_writeByName = 0;
	} else {
		try{
			_g2 = list.GetDoubleParameter(tmp_g2_name.get());
			_g2_writeByName = 1;
		} catch (BadParameter& ex){
			if(!_enable){
				_g2 = std::shared_ptr<DoubleParameter>(
						new DoubleParameter(
								tmp_g2_name.get(),0.0)
				);
				_g2_writeByName = 1;
			} else {
				BOOST_LOG_TRIVIAL(error) <<"AmpFlatteRes::Configure() | "
						"Requesting parameter "<<tmp_g2_name.get()<<" but"
						" was not found in parameter list. "
						"Quit since parameter is mandatory!";
				throw;
			}
		}
	}

	auto tmp_g2_massA = pt.get_optional<std::string>("g2_part1");
	if(!tmp_g2_massA)
		throw BadParameter("AmpFlatteRes::Configure() | "
				"g2_massA for "+_name+" not specified!");
	_g2_idA = tmp_g2_massA.get();
	_g2_massA = PhysConst::instance()->getMass(_g2_idA);

	auto tmp_g2_massB = pt.get_optional<std::string>("g2_part2");
	if(!tmp_g2_massB)
		throw BadParameter("AmpFlatteRes::Configure() | "
				"g2_massB for "+_name+" not specified!");
	_g2_idB = tmp_g2_massB.get();
	_g2_massB = PhysConst::instance()->getMass(_g2_idB);

	//Coupling to channel 3 (optional)
	auto tmp_g3_fix = pt.get<bool>("g3_fix",1);
	auto tmp_g3_min = pt.get<double>("g3_min",0.);
	auto tmp_g3_max = pt.get<double>("g3_max",10.);
	auto tmp_g3_name = pt.get_optional<std::string>("g3_name");
	if(!tmp_g3_name){
		auto tmp_g3= pt.get<double>("g3",0.0);
		_g3 = std::shared_ptr<DoubleParameter>(
				new DoubleParameter(
						"g3_"+_name,tmp_g3, tmp_g3_min, tmp_g3_max
				)
		);
		_g3->FixParameter(tmp_g3_fix);
		_g3_writeByName = 0;
		auto tmp_g3_massA = pt.get<std::string>("g3_part1","");
		_g3_idA = tmp_g3_massA;
		_g3_massA = PhysConst::instance()->getMass(tmp_g3_massA);
		auto tmp_g3_massB = pt.get<std::string>("g3_part2","");
		_g3_idB = tmp_g3_massB;
		_g3_massB = PhysConst::instance()->getMass(tmp_g3_massB);

		if(_enable && _g3->GetValue() != 0 ) {
		BOOST_LOG_TRIVIAL(info) << "AmpFlatteRes::Configure() | Adding 3rd channel: "
				<<_name<<" parName="<<_g3->GetName()<<" g3="<<_g3->GetValue()
				<<" g3_massA="<<_g3_massA<<" g3_massB="<<_g3_massB;
			list.AddParameter(_g3);
		}
	} else {
		try{
			_g3 = list.GetDoubleParameter(tmp_g3_name.get());
			_g3_writeByName = 1;
			auto tmp_g3_massA = pt.get<std::string>("g3_part1","");
			_g3_idA = tmp_g3_massA;
			_g3_massA = PhysConst::instance()->getMass(tmp_g3_massA);
			auto tmp_g3_massB = pt.get<std::string>("g3_part2","");
			_g3_idB = tmp_g3_massB;
			_g3_massB = PhysConst::instance()->getMass(tmp_g3_massB);
			BOOST_LOG_TRIVIAL(info) << "AmpFlatteRes::Configure() | Adding 3rd channel:"
					<<_name<<" parName="<<_g3->GetName()<<" g3="<<_g3->GetValue()
					<<" g3_massA="<<_g3_massA<<" g3_massB="<<_g3_massB;
		} catch (BadParameter& ex){
			if(!_enable){
				_g3 = std::shared_ptr<DoubleParameter>(
						new DoubleParameter(
								tmp_g3_name.get(),0.0)
				);
				_g3_writeByName = 1;
			} else {
				BOOST_LOG_TRIVIAL(error) <<"AmpFlatteRes::Configure() | "
						"Requesting parameter "<<tmp_g3_name.get()<<" but"
						" was not found in parameter list. "
						"Continue since parameter is not mandatory!";
				_g3 =std::shared_ptr<DoubleParameter>(
						new DoubleParameter("g3_"+_name,0.0)
				);
			}
		}
	}
	initialize();

	return;
}

void AmpFlatteRes::Save(boost::property_tree::ptree &pt)
{
	boost::property_tree::ptree amp;
	AmpAbsDynamicalFunction::put(amp);
	if(_g1_writeByName){
		amp.put("g1_name",_g1->GetName());
	} else {
		amp.put("g1", _g1->GetValue());
		amp.put("g1_fix", _g1->IsFixed());
		amp.put("g1_min", _g1->GetMinValue());
		amp.put("g1_max", _g1->GetMaxValue());
	}
	if(_g2_writeByName){
		amp.put("g2_name",_g2->GetName());
	} else {
		amp.put("g2", _g2->GetValue());
		amp.put("g2_fix", _g2->IsFixed());
		amp.put("g2_min", _g2->GetMinValue());
		amp.put("g2_max", _g2->GetMaxValue());
	}
	amp.put("g2_part1", _g2_idA);
	amp.put("g2_part2", _g2_idB);
	if(_g3->GetValue() != 0){
		if(_g3_writeByName){
			amp.put("g3_name",_g3->GetName());
		} else {
			amp.put("g3", _g3->GetValue());
			amp.put("g3_fix", _g3->IsFixed());
			amp.put("g3_min", _g3->GetMinValue());
			amp.put("g3_max", _g3->GetMaxValue());
		}
		amp.put("g3_part1", _g3_idA);
		amp.put("g3_part2", _g3_idB);
	}

	pt.add_child("Flatte", amp);

	return;
}

std::string AmpFlatteRes::to_str() const
{
	std::string dynAmp = AmpAbsDynamicalFunction::to_str();
	std::stringstream str;
	str<<_g1->to_str()<<std::endl;
	str<<_g2->to_str()<<std::endl;
	str<<_g3->to_str()<<std::endl;
	str<<"massB1="<<_g2_massA<<" massB2="<<_g2_massB;
	str<<" massC1="<<_g3_massA<<" massC2="<<_g3_massB<<std::endl;

	return dynAmp+str.str();
}

void AmpFlatteRes::CheckModified()
{
	AmpAbsDynamicalFunction::CheckModified();
	if( _g1->GetValue()!=tmp_g1 || _g2->GetValue()!=tmp_g2 || _g3->GetValue()!=tmp_g3 ) {
		SetModified();
		tmp_g1 = _g1->GetValue();
		tmp_g2 = _g2->GetValue();
		tmp_g3 = _g3->GetValue();
	}
	return;
}
std::complex<double> AmpFlatteRes::EvaluateAmp(dataPoint& point)
{
	CheckModified();

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double mSq=-999;
	switch(_subSys){
	case 3:
		mSq=(kin->getThirdVariableSq(point.getVal(0),point.getVal(1)));	break;
	case 4: mSq=(point.getVal(1)); break;
	case 5: mSq=(point.getVal(0)); break;
	}

	return dynamicalFunction(mSq,_mass->GetValue(),_ma,_mb,_g1->GetValue(),
			_g2_massA,_g2_massB,_g2->GetValue(),
			_g3_massA,_g3_massB,_g3->GetValue(),
			_spin,_mesonRadius->GetValue(), _ffType);
}

std::complex<double> AmpFlatteRes::dynamicalFunction(double mSq, double mR,
		double massA1, double massA2, double gA,
		double massB1, double massB2, double couplingB,
		double massC1, double massC2, double couplingC,
		unsigned int J, double mesonRadius, formFactorType ffType )
{
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
	if(couplingC!=0.0){
		//break-up momentum
		barrierC = Kinematics::FormFactor(sqrtS,massC1,massC2,J,mesonRadius, ffType)/
				Kinematics::FormFactor(mR,massC1,massC2,J,mesonRadius, ffType);
		gC = couplingC;
		//convert coupling to partial width of channel C
		gammaC = couplingToWidth(mSq,mR,gC,massC1,massC2,J,mesonRadius, ffType);
		//		qTermC = Kinematics::qValue(sqrtS,massC1,massC2) / Kinematics::qValue(mR,massC1,massC2);
		qTermC = std::complex<double>(1,0);
		termC = gammaC*barrierC*barrierC*std::pow(qTermC,(double)2*J+1);
	}

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
		std::cout<<"mpFlatteRes::dynamicalFunction() | "<<barrierA<<" "<<mR<<" "<<mSq<<" "
				<<massA1<<" "<<massA2<<std::endl;
		return 0;
	}
	return result;
}

std::shared_ptr<FunctionTree> AmpFlatteRes::SetupTree(
		allMasses& theMasses,allMasses& toyPhspSample,std::string suffix)
{
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double phspVol = kin->getPhspVolume();
	BOOST_LOG_TRIVIAL(info) << "AmpFlatteRes::setupBasicTree() | "
			<<_name << " nEvents=" <<theMasses.nEvents<<" nPhspEvents="
			<<toyPhspSample.nEvents;

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
	std::shared_ptr<AbsSquare> msqStrat(new AbsSquare(ParType::MDOUBLE));
	std::shared_ptr<MultAll> multDStrat(new MultAll(ParType::DOUBLE));
	std::shared_ptr<AddAll> addStrat(new AddAll(ParType::DOUBLE));
	std::shared_ptr<Complexify> complStrat(new Complexify(ParType::COMPLEX));
	std::shared_ptr<Inverse> invStrat(new Inverse(ParType::DOUBLE));
	std::shared_ptr<SquareRoot> sqRootStrat(new SquareRoot(ParType::DOUBLE));

	//----Add Nodes
	std::shared_ptr<FlatteStrategy> flatteStrat(
			new FlatteStrategy(_name,ParType::MCOMPLEX));
	std::shared_ptr<WignerDStrategy> angdStrat(
			new WignerDStrategy(_name,ParType::MDOUBLE) );

	newTree->createHead("Reso_"+_name, mmultStrat, theMasses.nEvents); //Reso=BW*C_*AD*N_
	newTree->createNode("C_"+_name, complStrat, "Reso_"+_name); //c
	newTree->createLeaf("Intens_"+_name, _mag, "C_"+_name); //r
	newTree->createLeaf("Phase_"+_name, _phase, "C_"+_name); //phi
	//Angular distribution
	newTree->insertTree(_wignerD.SetupTree(theMasses,suffix), "Reso_"+_name);

	//Flatte
	newTree->createNode("FlatteRes_"+_name, flatteStrat, "Reso_"+_name, theMasses.nEvents); //BW
	newTree->createLeaf("mass", _mass, "FlatteRes_"+_name);//use global parameter m0_a0 (asdfef)
	newTree->createLeaf("g1", _g1, "FlatteRes_"+_name);//use global parameter g1_a0 (asdfef)
	newTree->createLeaf("spin", _spin, "FlatteRes_"+_name); //spin
	newTree->createLeaf("mesonRadius", _mesonRadius , "FlatteRes_"+_name); //d
	newTree->createLeaf("formFactorType", _ffType , "FlatteRes_"+_name); //d
	newTree->createLeaf("subSysFlag", _subSys, "FlatteRes_"+_name); //_subSysFlag
	switch(_subSys){
	case 3:{ //reso in sys of particles 1&2
		newTree->createLeaf("ma_"+_name, kin->m1, "FlatteRes_"+_name); //ma
		newTree->createLeaf("mb_"+_name, kin->m2, "FlatteRes_"+_name); //mb
		break;
	} case 4: { //reso in sys of particles 1&3
		newTree->createLeaf("ma_"+_name, kin->m1, "FlatteRes_"+_name); //ma
		newTree->createLeaf("mb_"+_name, kin->m3, "FlatteRes_"+_name); //mb
		break;
	} case 5: { //reso in sys of particles 2&3
		newTree->createLeaf("ma_"+_name, kin->m2, "FlatteRes_"+_name); //ma
		newTree->createLeaf("mb_"+_name, kin->m3, "FlatteRes_"+_name); //mb
		break;
	} default: {
		BOOST_LOG_TRIVIAL(error)<<"AmpFlatteRes::setupBasicTree() | "
				"Subsys not found!";
	}
	}
	newTree->createLeaf("g2", _g2, "FlatteRes_"+_name);
	newTree->createLeaf("massB1", _g2_massA, "FlatteRes_"+_name);
	newTree->createLeaf("massB2", _g2_massB, "FlatteRes_"+_name);
	newTree->createLeaf("g3", _g3, "FlatteRes_"+_name);
	newTree->createLeaf("massC1", _g3_massA, "FlatteRes_"+_name);
	newTree->createLeaf("massC2", _g3_massB, "FlatteRes_"+_name);
	newTree->createLeaf("m12sq", m12sq, "FlatteRes_"+_name); //mc
	newTree->createLeaf("m13sq", m13sq, "FlatteRes_"+_name); //mb
	newTree->createLeaf("m23sq", m23sq, "FlatteRes_"+_name); //ma

	//Normalization
	if(_normStyle==normStyle::none){
		newTree->createLeaf("N_"+_name, 1., "Reso_"+_name);
	} else {
		newTree->createNode("N_"+_name, sqRootStrat, "Reso_"+_name); //N = sqrt(NSq)
		newTree->createNode("NSq_"+_name, multDStrat, "N_"+_name); //NSq = N_phspMC * 1/PhspVolume * 1/Sum(|A|^2)
		newTree->createLeaf("PhspSize_"+_name, toyPhspSample.nEvents, "NSq_"+_name); // N_phspMC
		newTree->createLeaf("PhspVolume_"+_name, 1/phspVol, "NSq_"+_name); // 1/PhspVolume
		newTree->createNode("InvSum_"+_name, invStrat, "NSq_"+_name); //1/Sum(|A|^2)
		newTree->createNode("Sum_"+_name, addStrat, "InvSum_"+_name); //Sum(|A|^2)
		newTree->createNode("AbsVal_"+_name, msqStrat, "Sum_"+_name); //|A_i|^2

		newTree->createNode("NormReso_"+_name, mmultStrat, "AbsVal_"+_name,
				toyPhspSample.nEvents); //BW
		//Angular distribution (Normalization)
		newTree->insertTree(_wignerD.SetupTree(toyPhspSample,suffix),
				"NormReso_"+_name);
		//Flatte (Normalization)
		newTree->createNode("NormFlatte_"+_name, flatteStrat, "NormReso_"+_name,
				toyPhspSample.nEvents); //BW
		newTree->createLeaf("mass", _mass, "NormFlatte_"+_name);//use local parameter m0_a0
		newTree->createLeaf("g1", _g1, "NormFlatte_"+_name);//use global parameter g1_a0 (asdfef)
		newTree->createLeaf("spin", _spin, "NormFlatte_"+_name); //spin
		newTree->createLeaf("mesonRadius",  _mesonRadius, "NormFlatte_"+_name); //d
		newTree->createLeaf("formFactorType", _ffType , "NormFlatte_"+_name); //d
		newTree->createLeaf("subSysFlag", _subSys, "NormFlatte_"+_name); //_subSysFlag
		switch(_subSys){
		case 3:{ //reso in sys of particles 1&2
			newTree->createLeaf("ma_"+_name, kin->m1, "NormFlatte_"+_name); //ma
			newTree->createLeaf("mb_"+_name, kin->m2, "NormFlatte_"+_name); //mb
			break;
		} case 4:{ //reso in sys of particles 1&3
			newTree->createLeaf("ma_"+_name, kin->m1, "NormFlatte_"+_name); //ma
			newTree->createLeaf("mb_"+_name, kin->m3, "NormFlatte_"+_name); //mb
			break;
		} case 5:{ //reso in sys of particles 2&3
			newTree->createLeaf("ma_"+_name, kin->m2, "NormFlatte_"+_name); //ma
			newTree->createLeaf("mb_"+_name, kin->m3, "NormFlatte_"+_name); //mb
			break;
		}
		default:{
			BOOST_LOG_TRIVIAL(error)<<"AmpSumIntensity::setupBasicTree(): Subsys not found!!";
		}
		}
		newTree->createLeaf("g2", _g2, "NormFlatte_"+_name);
		newTree->createLeaf("massB1", _g2_massA, "NormFlatte_"+_name);
		newTree->createLeaf("massB2", _g2_massB, "NormFlatte_"+_name);
		newTree->createLeaf("g3", _g3, "NormFlatte_"+_name);
		newTree->createLeaf("massC1", _g3_massA, "NormFlatte_"+_name);
		newTree->createLeaf("massC2", _g3_massB, "NormFlatte_"+_name);
		newTree->createLeaf("m12sq_phsp", m12sq_phsp, "NormFlatte_"+_name);
		newTree->createLeaf("m13sq_phsp", m13sq_phsp, "NormFlatte_"+_name);
		newTree->createLeaf("m23sq_phsp", m23sq_phsp, "NormFlatte_"+_name);
	}

	return newTree;
}

bool FlatteStrategy::execute(ParameterList& paras,
		std::shared_ptr<AbsParameter>& out)
{
	//Debug
	//	BOOST_LOG_TRIVIAL(debug) <<"FlatteStrategy::execute() | start";
	//	for( int i=0; i<paras.GetDoubleParameters().size(); ++i){
	//		BOOST_LOG_TRIVIAL(debug)<<paras.GetDoubleParameter(i)->GetName();
	//	}

	//Check parameter type
	if( checkType != out->type() )
		throw( WrongParType(	std::string("Output Type ")
	+ParNames[out->type()] + std::string(" conflicts expected type ")
	+ParNames[checkType]+std::string(" of ")+name+" BW strat")
		);

	//Check size of parameter list
	if( paras.GetNDouble() != 14 && paras.GetNDouble() != 17)
		throw( BadParameter("FlatteStrategy::execute() | "
				"number of DoubleParameters does not match!")
		);

	double m0, d, ma, mb, g1, g2, massB1, massB2, g3, massC1, massC2;
	unsigned int spin, subSys;
	int ffType;
	/** Get parameters from ParameterList:
	 * We use the same order of the parameters as was used during tree
	 * construction
	 */
	m0 = paras.GetDoubleParameter(0)->GetValue();
	g1 = paras.GetDoubleParameter(1)->GetValue();
	spin = (unsigned int)paras.GetDoubleParameter(2)->GetValue();
	d = paras.GetDoubleParameter(3)->GetValue();
	ffType = paras.GetDoubleParameter(4)->GetValue();
	subSys = paras.GetDoubleParameter(5)->GetValue();
	ma = paras.GetDoubleParameter(6)->GetValue();
	mb = paras.GetDoubleParameter(7)->GetValue();
	g2 = paras.GetDoubleParameter(8)->GetValue();
	massB1 = paras.GetDoubleParameter(9)->GetValue();
	massB2 = paras.GetDoubleParameter(10)->GetValue();
	g3 = paras.GetDoubleParameter(11)->GetValue();
	massC1 = paras.GetDoubleParameter(12)->GetValue();
	massC2 = paras.GetDoubleParameter(13)->GetValue();

	//	BOOST_LOG_TRIVIAL(debug) << "FlatteStrategy::execute() | mR="<<m0
	//			<<" g1="<<g1<<" spin="<<spin<<" radius="<<d<<" ffType="<<ffType
	//			<<" subSys="<<subSys<<" mA1="<<ma<<" mA2="<<mb
	//			<<" g2="<<g2<<" mB1="<<massB1<<" mB2="<<massB2;
	//			<<" g3="<<g3<<" mC1="<<massC1<<" mC2="<<massC2;

	//MultiDim output, must have multidim Paras in input
	if(checkType == ParType::MCOMPLEX){
		std::shared_ptr<MultiDouble> mp;
		try {
			switch(subSys){
			case 3:{ mp  = (paras.GetMultiDouble(0)); break; }
			case 4:{ mp  = (paras.GetMultiDouble(1)); break; }
			case 5:{ mp  = (paras.GetMultiDouble(2)); break; }
			}
		} catch (std::exception &ex) {
			BOOST_LOG_TRIVIAL(error) << "FlatteStrategy::execute() | "
					<<ex.what();
			throw(WrongParType("FlatteStrategy::execute() | "
					"Failed to obtain data vector from parameter list!")
			);
		}
		std::vector<std::complex<double> > results(mp->GetNValues(),
				std::complex<double>(0.));
		//calc BW for each point
		for(unsigned int ele=0; ele<mp->GetNValues(); ele++){
			double mSq = (mp->GetValue(ele));
			try{
				results[ele] = AmpFlatteRes::dynamicalFunction(
						mSq,m0,ma,mb,g1,massB1,massB2,g2,massC1,massC2,g3,
						spin,d, formFactorType(ffType)
				);
			} catch (std::exception& ex) {
				BOOST_LOG_TRIVIAL(error) << "FlatteStrategy::execute() | "
						<<ex.what();
				throw( std::runtime_error("FlatteStrategy::execute() | "
						"Evaluation of dynamic function failed!")
				);
			}
		}
		out = std::shared_ptr<AbsParameter>(
				new MultiComplex(out->GetName(),results));
		//		BOOST_LOG_TRIVIAL(debug) <<"FlatteStrategy::execute() | finished!";
		return true;
	}//end multicomplex output

	double mSq;
	try {
		switch(subSys){
		case 3:{ mSq  = paras.GetDoubleParameter(14)->GetValue(); break; }
		case 4:{ mSq  = paras.GetDoubleParameter(15)->GetValue(); break; }
		case 5:{ mSq  = paras.GetDoubleParameter(16)->GetValue(); break; }
		}
	} catch (std::exception &ex) {
		BOOST_LOG_TRIVIAL(error) << "FlatteStrategy::execute() | "
				<<ex.what();
		throw(WrongParType("FlatteStrategy::execute() | "
				"Failed to obtain data from parameter list!")
		);
	}
	std::complex<double> result;
	try{
		result = AmpFlatteRes::dynamicalFunction(
				mSq,m0,ma,mb,g1,massB1,massB2,g2,massC1,massC2,g3,
				spin,d, formFactorType(ffType)
		);
	} catch (std::exception& ex) {
		BOOST_LOG_TRIVIAL(error) << "FlatteStrategy::execute() | "
				<<ex.what();
		throw( std::runtime_error("FlatteStrategy::execute() | "
				"Evaluation of dynamic function failed!")
		);
	}
	out = std::shared_ptr<AbsParameter>(
			new ComplexParameter(out->GetName(), result));
	//	BOOST_LOG_TRIVIAL(debug) <<"FlatteStrategy::execute() | finished!";
	return true;
}




