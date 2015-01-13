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

#include <vector>
#include <memory>
#include <ctime>
#include <exception>
#include <stdexcept>
#include <stdlib.h>

#include "Core/PhysConst.hpp"
#include "Core/Functions.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Physics/AmplitudeSum/AmpSumOfAmplitudes.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"

#include "gsl/gsl_monte.h"
#include "gsl/gsl_monte_plain.h"
#include "gsl/gsl_monte_miser.h"
#include "gsl/gsl_monte_vegas.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include "boost/function.hpp"
#include <boost/log/trivial.hpp>
using namespace boost::log;


AmpSumIntensity::AmpSumIntensity(const AmpSumIntensity& other) : nAmps(other.nAmps), _dpArea(other._dpArea),
		_normStyle(other._normStyle),maxVal(other.maxVal),ampSetup(other.ampSetup),
		totAmp(other.totAmp), _calcMaxFcnVal(other._calcMaxFcnVal), _maxFcnVal(other._maxFcnVal),
		_nCalls(other._nCalls){
}

AmpSumIntensity::AmpSumIntensity(AmplitudeSetup ini, normStyle ns, std::shared_ptr<Efficiency> eff,
		unsigned int nCalls) :
							totAmp("relBWsumAmplitude"), ampSetup(ini),
							_normStyle(ns), _calcNorm(1), _dpArea(1.),
							_calcMaxFcnVal(0),eff_(eff),_nCalls(nCalls)
{
	init();
}

AmpSumIntensity::AmpSumIntensity(AmplitudeSetup ini, std::shared_ptr<Efficiency> eff,
		unsigned int nCalls) :
							totAmp("relBWsumAmplitude"), ampSetup(ini),
							_normStyle(none), _calcNorm(0), _dpArea(1.),
							_calcMaxFcnVal(0),eff_(eff),_nCalls(nCalls)
{
	init();
}

void AmpSumIntensity::init(){
	result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("AmpSumResult")));

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	_calcNorm=1;
	if(_normStyle==normStyle::none) _calcNorm=0;

	params.AddParameter( std::shared_ptr<DoubleParameter> (
			new DoubleParameter("motherRadius",1.5) ));
	params.GetDoubleParameter("motherRadius")->FixParameter(1);
	/* For peter's analysis the a_0+ and a_00 share the same coupling. To implement
	 * this in the model, we have to do the following work-a-round.
	 * Don't forget to adjust the iterator 'paramsPos' and 'g1Itr', if you comment out that line.
	 * Search for "(asdfef)" to find all positions that have to be adjusted
	 */
	params.AddParameter( std::shared_ptr<DoubleParameter> (
			new DoubleParameter("g1_a_0",0.464) )); //(asdfef)

	for(std::vector<Resonance>::iterator reso=ampSetup.getResonances().begin(); reso!=ampSetup.getResonances().end(); reso++){
		Resonance tmp = (*reso);
		if(!tmp.m_enable) continue;
		std::string name = tmp.m_name;
		namer.push_back(name);
		unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("mag_"+name,tmp.m_strength,tmp.m_strength_min,tmp.m_strength_max) ) );
		params.GetDoubleParameter("mag_"+name)->FixParameter(tmp.m_strength_fix);
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("phase_"+name,tmp.m_phase,tmp.m_phase_min,tmp.m_phase_max) ) );
		params.GetDoubleParameter("phase_"+name)->FixParameter(tmp.m_phase_fix);
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("m0_"+name,tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max) ) );
		params.GetDoubleParameter("m0_"+name)->FixParameter(tmp.m_mass_fix);
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("width_"+name,tmp.m_width, tmp.m_width_min, tmp.m_width_max) ) );
		params.GetDoubleParameter("width_"+name)->FixParameter(tmp.m_width_fix);
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("d_"+name,tmp.m_mesonRadius) ) );
		params.GetDoubleParameter("d_"+name)->FixParameter(1);

		std::shared_ptr<AmpRelBreitWignerRes> tmpbw(new AmpRelBreitWignerRes(name.c_str(),
				params.GetDoubleParameter("m0_"+name), params.GetDoubleParameter("width_"+name),
				params.GetDoubleParameter("d_"+name),
				params.GetDoubleParameter("motherRadius"),
				subSys, tmp.m_spin,tmp.m_m,tmp.m_n) );
		totAmp.addBW(tmpbw, params.GetDoubleParameter("mag_"+name), params.GetDoubleParameter("phase_"+name));

		//setting normalization between amplitudes
		double norm=1.0;
		if(norm<0 || _calcNorm) {//recalculate normalization
			norm = normReso(tmpbw);
		}
		tmpbw->SetNormalization(1/norm);
	}// end loop over resonances

	for(std::vector<ResonanceFlatte>::iterator reso=ampSetup.getResonancesFlatte().begin(); reso!=ampSetup.getResonancesFlatte().end(); reso++){
		ResonanceFlatte tmp = (*reso);
		if(!tmp.m_enable) continue;
		std::string name = tmp.m_name;
		namer.push_back(name);
		unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("mag_"+name,tmp.m_strength,tmp.m_strength_min,tmp.m_strength_max) ) );
		params.GetDoubleParameter("mag_"+name)->FixParameter(tmp.m_strength_fix);
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("phase_"+name,tmp.m_phase,tmp.m_phase_min,tmp.m_phase_max) ) );
		params.GetDoubleParameter("phase_"+name)->FixParameter(tmp.m_phase_fix);
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("m0_"+name,tmp.m_mass, tmp.m_mass_min, tmp.m_mass_max)  ));
		params.GetDoubleParameter("m0_"+name)->FixParameter(tmp.m_mass_fix);
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("g1_"+name,tmp.m_g1,tmp.m_g1_min,tmp.m_g1_max) ) );
		params.GetDoubleParameter("g1_"+name)->FixParameter(tmp.m_g1_fix);
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("g2_"+name,tmp.m_g2) ) );
		params.GetDoubleParameter("g2_"+name)->FixParameter(1);
		params.AddParameter( std::shared_ptr<DoubleParameter> (
				new DoubleParameter("d_"+name,tmp.m_mesonRadius) ) );
		params.GetDoubleParameter("d_"+name)->FixParameter(1);

		std::shared_ptr<AmpFlatteRes> tmpbw(new AmpFlatteRes(name.c_str(),
				params.GetDoubleParameter("m0_"+name), params.GetDoubleParameter("d_"+name),
				params.GetDoubleParameter("motherRadius"),
//				params.GetDoubleParameter("g1_"+name),params.GetDoubleParameter("g2_"+name),
				params.GetDoubleParameter("g1_a_0"),params.GetDoubleParameter("g2_"+name),
				PhysConst::instance()->getMass(tmp.m_g2_part1),
				PhysConst::instance()->getMass(tmp.m_g2_part2),
				subSys, tmp.m_spin,tmp.m_m,tmp.m_n) );

		totAmp.addBW(tmpbw, params.GetDoubleParameter("mag_"+name), params.GetDoubleParameter("phase_"+name));

		double norm=1.0;
		if(norm<0 || _calcNorm)	norm = normReso(tmpbw);
		tmpbw->SetNormalization(1./norm);
	}// end loop over resonancesFlatte

	nAmps=namer.size();
	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity: completed setup!";
}

std::shared_ptr<FunctionTree> AmpSumIntensity::setupBasicTree(
		allMasses& theMasses,allMasses& toyPhspSample,std::string suffix)
{
	BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::setupTree() generating new tree!";
	if(theMasses.nEvents==0){
		BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::setupBasicTree() data sample empty!";
		return std::shared_ptr<FunctionTree>();
	}
	if(toyPhspSample.nEvents==0){
		BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::setupBasicTree() toy sample empty!";
		return std::shared_ptr<FunctionTree>();
	}
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	_dpArea = kin->getPhspVolume();

	//------------Setup Tree---------------------
	std::shared_ptr<FunctionTree> newTree = std::shared_ptr<FunctionTree>(new FunctionTree());
	//------------Setup Tree Pars---------------------
	std::shared_ptr<MultiDouble> m23sq = std::shared_ptr<MultiDouble>( new MultiDouble("m23sq",theMasses.masses_sq.at( std::make_pair(2,3) )) );
	std::shared_ptr<MultiDouble> m13sq = std::shared_ptr<MultiDouble>( new MultiDouble("m13sq",theMasses.masses_sq.at( std::make_pair(1,3) )) );
	std::shared_ptr<MultiDouble> m12sq = std::shared_ptr<MultiDouble>( new MultiDouble("m12sq",theMasses.masses_sq.at( std::make_pair(1,2) )) );
	std::shared_ptr<MultiDouble> m23sq_phsp = std::shared_ptr<MultiDouble>( new MultiDouble("m23sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(2,3) )) );
	std::shared_ptr<MultiDouble> m13sq_phsp = std::shared_ptr<MultiDouble>( new MultiDouble("m13sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(1,3) )) );
	std::shared_ptr<MultiDouble> m12sq_phsp = std::shared_ptr<MultiDouble>( new MultiDouble("m12sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(1,2) )) );

	//----Strategies needed
	std::shared_ptr<MultAll> mmultStrat = std::shared_ptr<MultAll>(new MultAll(ParType::MCOMPLEX));
	std::shared_ptr<MultAll> mmultDStrat = std::shared_ptr<MultAll>(new MultAll(ParType::MDOUBLE));
	std::shared_ptr<AddAll> maddStrat = std::shared_ptr<AddAll>(new AddAll(ParType::MCOMPLEX));
	std::shared_ptr<AbsSquare> msqStrat = std::shared_ptr<AbsSquare>(new AbsSquare(ParType::MDOUBLE));
	std::shared_ptr<LogOf> mlogStrat = std::shared_ptr<LogOf>(new LogOf(ParType::MDOUBLE));
	std::shared_ptr<MultAll> multStrat = std::shared_ptr<MultAll>(new MultAll(ParType::COMPLEX));
	std::shared_ptr<MultAll> multDStrat = std::shared_ptr<MultAll>(new MultAll(ParType::DOUBLE));
	std::shared_ptr<AddAll> addStrat = std::shared_ptr<AddAll>(new AddAll(ParType::DOUBLE));
	std::shared_ptr<AddAll> addComplexStrat = std::shared_ptr<AddAll>(new AddAll(ParType::COMPLEX));
	std::shared_ptr<AbsSquare> sqStrat = std::shared_ptr<AbsSquare>(new AbsSquare(ParType::DOUBLE));
	std::shared_ptr<LogOf> logStrat = std::shared_ptr<LogOf>(new LogOf(ParType::DOUBLE));
	std::shared_ptr<Complexify> complStrat = std::shared_ptr<Complexify>(new Complexify(ParType::COMPLEX));
	std::shared_ptr<Inverse> invStrat = std::shared_ptr<Inverse>(new Inverse(ParType::DOUBLE));
	std::shared_ptr<SquareRoot> sqRootStrat = std::shared_ptr<SquareRoot>(new SquareRoot(ParType::DOUBLE));

	newTree->createHead("Amplitude"+suffix, maddStrat, theMasses.nEvents);

	//----Add Resonances
	for(std::vector<Resonance>::iterator reso=ampSetup.getResonances().begin(); reso!=ampSetup.getResonances().end(); reso++){
		Resonance tmp = (*reso);
		std::string name = tmp.m_name;
		if(!tmp.m_enable) continue;
		BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::setupTree() adding "<<tmp.m_name<<" (BW) to tree.";

		//----Add Nodes
		std::shared_ptr<BreitWignerStrategy> rbwStrat =
				std::shared_ptr<BreitWignerStrategy>(
						new BreitWignerStrategy(tmp.m_name,ParType::MCOMPLEX) );
		std::shared_ptr<BreitWignerPhspStrategy> rbwPhspStrat =
				std::shared_ptr<BreitWignerPhspStrategy>(
						new BreitWignerPhspStrategy(tmp.m_name,ParType::MCOMPLEX) );
		std::shared_ptr<WignerDStrategy> angdStrat =
				std::shared_ptr<WignerDStrategy>(
						new WignerDStrategy(tmp.m_name,ParType::MDOUBLE) );

		unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;
		newTree->createNode("Reso_"+tmp.m_name, mmultStrat, "Amplitude"+suffix, theMasses.nEvents); //Reso=BW*C_*AD*N_
		newTree->createNode("BW_"+tmp.m_name, mmultStrat , "Reso_"+tmp.m_name, theMasses.nEvents); //BW
		newTree->createNode("RelBW_"+tmp.m_name, rbwStrat, "BW_"+tmp.m_name, theMasses.nEvents); //BW
		newTree->createNode("C_"+tmp.m_name, complStrat, "Reso_"+tmp.m_name); //m0c
		newTree->createLeaf("Intens_"+tmp.m_name, params.GetDoubleParameter("mag_"+name), "C_"+tmp.m_name); //r
		newTree->createLeaf("Phase_"+tmp.m_name, params.GetDoubleParameter("phase_"+name), "C_"+tmp.m_name); //phi
		newTree->createNode("AngD_"+tmp.m_name, angdStrat, "Reso_"+tmp.m_name, theMasses.nEvents); //AD

		//Breit-Wigner
		newTree->createLeaf("m0_"+tmp.m_name, params.GetDoubleParameter("m0_"+name), "RelBW_"+tmp.m_name); //m0
		newTree->createLeaf("m23sq", m23sq, "RelBW_"+tmp.m_name); //ma
		newTree->createLeaf("m13sq", m13sq, "RelBW_"+tmp.m_name); //mb
		newTree->createLeaf("m12sq", m12sq, "RelBW_"+tmp.m_name); //mc
		newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "RelBW_"+tmp.m_name); //subSysFlag
		newTree->createLeaf("spin_"+tmp.m_name, tmp.m_spin, "RelBW_"+tmp.m_name); //spin
		newTree->createLeaf("d_"+tmp.m_name, params.GetDoubleParameter("d_"+name), "RelBW_"+tmp.m_name); //d
		newTree->createLeaf("width_"+tmp.m_name, params.GetDoubleParameter("width_"+name), "RelBW_"+tmp.m_name); //resWidth
		//Angular distribution
		newTree->createLeaf("m23sq", m23sq, "AngD_"+tmp.m_name); //ma
		newTree->createLeaf("m13sq", m13sq, "AngD_"+tmp.m_name); //mb
		newTree->createLeaf("m12sq", m12sq, "AngD_"+tmp.m_name); //mc
		newTree->createLeaf("M", kin->M, "AngD_"+tmp.m_name); //M
		newTree->createLeaf("m1", kin->m1, "AngD_"+tmp.m_name); //m1
		newTree->createLeaf("m2", kin->m2, "AngD_"+tmp.m_name); //m2
		newTree->createLeaf("m3", kin->m3, "AngD_"+tmp.m_name); //m3
		newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "AngD_"+tmp.m_name); //subSysFlag
		newTree->createLeaf("spin_"+tmp.m_name,tmp.m_spin, "AngD_"+tmp.m_name); //spin
		newTree->createLeaf("m_"+tmp.m_name, tmp.m_m, "AngD_"+tmp.m_name); //OutSpin 1
		newTree->createLeaf("n_"+tmp.m_name, tmp.m_n, "AngD_"+tmp.m_name); //OutSpin 2

		//adding nodes and leafs for calculation of normalization
		if(_normStyle==normStyle::none){
			newTree->createLeaf("N_"+tmp.m_name, 1., "BW_"+tmp.m_name);
		}else{
			//Normalization parameter for dynamical amplitude
			newTree->createNode("N_"+tmp.m_name, sqRootStrat, "BW_"+tmp.m_name); //N = sqrt(NSq)
			newTree->createNode("NSq_"+tmp.m_name, multDStrat, "N_"+tmp.m_name); //NSq = (PhspVolume/N_phspMC * Sum(|A|^2))^-1
			newTree->createLeaf("PhspSize_"+tmp.m_name, toyPhspSample.nEvents, "NSq_"+tmp.m_name); // N_phspMC
			newTree->createLeaf("PhspVolume_"+tmp.m_name, 1/_dpArea, "NSq_"+tmp.m_name); // 1/PhspVolume
			newTree->createNode("InvSum_"+tmp.m_name, invStrat, "NSq_"+tmp.m_name); //1/Sum(|A|^2)
			newTree->createNode("Sum_"+tmp.m_name, addStrat, "InvSum_"+tmp.m_name); //Sum(|A|^2)
			newTree->createNode("AbsVal_"+tmp.m_name, msqStrat, "Sum_"+tmp.m_name); //|A_i|^2
			//Breit-Wigner (Normalization)
			newTree->createNode("NormBW_"+tmp.m_name, rbwPhspStrat, "AbsVal_"+tmp.m_name, toyPhspSample.nEvents); //BW
			newTree->createLeaf("m0_"+tmp.m_name, params.GetDoubleParameter("m0_"+name), "NormBW_"+tmp.m_name); //m0
			newTree->createLeaf("m23sq_phsp", m23sq_phsp, "NormBW_"+tmp.m_name); //ma
			newTree->createLeaf("m13sq_phsp", m13sq_phsp, "NormBW_"+tmp.m_name); //mb
			newTree->createLeaf("m12sq_phsp", m12sq_phsp, "NormBW_"+tmp.m_name); //mc
			newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "NormBW_"+tmp.m_name); //subSysFlag
			newTree->createLeaf("spin_"+tmp.m_name, tmp.m_spin, "NormBW_"+tmp.m_name); //spin
			newTree->createLeaf("d_"+tmp.m_name, params.GetDoubleParameter("d_"+name), "NormBW_"+tmp.m_name); //d
			newTree->createLeaf("width_"+tmp.m_name, params.GetDoubleParameter("width_"+name), "NormBW_"+tmp.m_name); //resWidth
		}
		switch(subSys){
		case 3:{ //reso in sys of particles 1&2
			newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "RelBW_"+tmp.m_name); //ma
			newTree->createLeaf("mb_"+tmp.m_name, kin->m2, "RelBW_"+tmp.m_name); //mb
			if(_normStyle!=normStyle::none){
				newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "NormBW_"+tmp.m_name); //ma
				newTree->createLeaf("mb_"+tmp.m_name, kin->m2, "NormBW_"+tmp.m_name); //mb
			}
			break;
		}
		case 4:{ //reso in sys of particles 1&3
			//newTree->createLeaf("mym_"+tmp.m_name, m13, "RelBW_"+tmp.m_name); //m
			newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "RelBW_"+tmp.m_name); //ma
			newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "RelBW_"+tmp.m_name); //mb
			if(_normStyle!=normStyle::none){
				newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "NormBW_"+tmp.m_name); //ma
				newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "NormBW_"+tmp.m_name); //mb
			}
			break;
		}
		case 5:{ //reso in sys of particles 2&3
			//newTree->createLeaf("mym_"+tmp.m_name, m23, "RelBW_"+tmp.m_name); //m
			newTree->createLeaf("ma_"+tmp.m_name, kin->m2, "RelBW_"+tmp.m_name); //ma
			newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "RelBW_"+tmp.m_name); //mb
			if(_normStyle!=normStyle::none){
				newTree->createLeaf("ma_"+tmp.m_name, kin->m2, "NormBW_"+tmp.m_name); //ma
				newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "NormBW_"+tmp.m_name); //mb
			}
			break;
		}
		default:{
			BOOST_LOG_TRIVIAL(error)<<"AmpSumIntensity::setupBasicTree(): Subsys not found!!";
		}
		}
	}// end loop over resonances

	for(std::vector<ResonanceFlatte>::iterator reso=ampSetup.getResonancesFlatte().begin(); reso!=ampSetup.getResonancesFlatte().end(); reso++){
		ResonanceFlatte tmp = (*reso);
		std::string name = tmp.m_name;
		if(!tmp.m_enable) continue;
		BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::setupBasicTree() adding "<<tmp.m_name<<" (Flatte) to tree.";

		//----Add Nodes
		std::shared_ptr<FlatteStrategy> flatteStrat = std::shared_ptr<FlatteStrategy>(new FlatteStrategy(tmp.m_name,ParType::MCOMPLEX));
		std::shared_ptr<FlattePhspStrategy> flattePhspStrat = std::shared_ptr<FlattePhspStrategy>(new FlattePhspStrategy(tmp.m_name,ParType::MCOMPLEX));
		std::shared_ptr<WignerDStrategy> angdStrat = std::shared_ptr<WignerDStrategy>(new WignerDStrategy(tmp.m_name,ParType::MDOUBLE));
		//		std::shared_ptr<WignerDphspStrategy> angdPhspStrat = std::shared_ptr<WignerDphspStrategy>(new WignerDphspStrategy(tmp.m_name,ParType::MDOUBLE));
		unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;
		newTree->createNode("Reso_"+tmp.m_name, mmultStrat, "Amplitude"+suffix, theMasses.nEvents); //Reso=BW*C_*AD*N_
		newTree->createNode("Flatte_"+tmp.m_name, mmultStrat , "Reso_"+tmp.m_name, theMasses.nEvents); //BW
		newTree->createNode("FlatteRes_"+tmp.m_name, flatteStrat, "Flatte_"+tmp.m_name, theMasses.nEvents); //BW
		newTree->createNode("C_"+tmp.m_name, complStrat, "Reso_"+tmp.m_name); //c
		newTree->createLeaf("Intens_"+tmp.m_name, params.GetDoubleParameter("mag_"+name), "C_"+tmp.m_name); //r
		newTree->createLeaf("Phase_"+tmp.m_name, params.GetDoubleParameter("phase_"+name), "C_"+tmp.m_name); //phi
		newTree->createNode("AngD_"+tmp.m_name, angdStrat, "Reso_"+tmp.m_name, theMasses.nEvents); //AD

		//Flatte
		newTree->createLeaf("m0_"+tmp.m_name, params.GetDoubleParameter("m0_"+name), "FlatteRes_"+tmp.m_name); //m0
		newTree->createLeaf("m23sq", m23sq, "FlatteRes_"+tmp.m_name); //ma
		newTree->createLeaf("m13sq", m13sq, "FlatteRes_"+tmp.m_name); //mb
		newTree->createLeaf("m12sq", m12sq, "FlatteRes_"+tmp.m_name); //mc
		newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "FlatteRes_"+tmp.m_name); //subSysFlag
		newTree->createLeaf("spin_"+tmp.m_name, tmp.m_spin, "FlatteRes_"+tmp.m_name); //spin
		newTree->createLeaf("d_"+tmp.m_name, params.GetDoubleParameter("d_"+name) , "FlatteRes_"+tmp.m_name); //d
		newTree->createLeaf("mHiddenA_"+tmp.m_name, \
				PhysConst::instance()->getMass(tmp.m_g2_part1), "FlatteRes_"+tmp.m_name);
		newTree->createLeaf("mHiddenB_"+tmp.m_name, \
				PhysConst::instance()->getMass(tmp.m_g2_part2), "FlatteRes_"+tmp.m_name);
//		newTree->createLeaf("g1_a_0", params.GetDoubleParameter("g1_a_0"), "FlatteRes_"+tmp.m_name);//use global parameter g1_a0 (asdfef)
		newTree->createLeaf("g1_a_0", params.GetDoubleParameter("g1_"+name), "FlatteRes_"+tmp.m_name);//use local parameter g1_a0
		newTree->createLeaf("g2_"+tmp.m_name, params.GetDoubleParameter("g2_"+name), "FlatteRes_"+tmp.m_name);
		//Angular distribution
		newTree->createLeaf("m23sq", m23sq, "AngD_"+tmp.m_name); //ma
		newTree->createLeaf("m13sq", m13sq, "AngD_"+tmp.m_name); //mb
		newTree->createLeaf("m12sq", m12sq, "AngD_"+tmp.m_name); //mc
		newTree->createLeaf("M", kin->M, "AngD_"+tmp.m_name); //M
		newTree->createLeaf("m1", kin->m1, "AngD_"+tmp.m_name); //m1
		newTree->createLeaf("m2", kin->m2, "AngD_"+tmp.m_name); //m2
		newTree->createLeaf("m3", kin->m3, "AngD_"+tmp.m_name); //m3
		newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "AngD_"+tmp.m_name); //subSysFlag
		newTree->createLeaf("spin_"+tmp.m_name,tmp.m_spin, "AngD_"+tmp.m_name); //spin
		newTree->createLeaf("m_"+tmp.m_name, tmp.m_m, "AngD_"+tmp.m_name); //OutSpin 1
		newTree->createLeaf("n_"+tmp.m_name, tmp.m_n, "AngD_"+tmp.m_name); //OutSpin 2

		//Normalization
		if(_normStyle!=normStyle::none){
			newTree->createNode("N_"+tmp.m_name, sqRootStrat, "Flatte_"+tmp.m_name); //N = sqrt(NSq)
			newTree->createNode("NSq_"+tmp.m_name, multDStrat, "N_"+tmp.m_name); //NSq = N_phspMC * 1/PhspVolume * 1/Sum(|A|^2)
			newTree->createLeaf("PhspSize_"+tmp.m_name, toyPhspSample.nEvents, "NSq_"+tmp.m_name); // N_phspMC
			newTree->createLeaf("PhspVolume_"+tmp.m_name, 1/_dpArea, "NSq_"+tmp.m_name); // 1/PhspVolume
			newTree->createNode("InvSum_"+tmp.m_name, invStrat, "NSq_"+tmp.m_name); //1/Sum(|A|^2)
			newTree->createNode("Sum_"+tmp.m_name, addStrat, "InvSum_"+tmp.m_name); //Sum(|A|^2)
			newTree->createNode("AbsVal_"+tmp.m_name, msqStrat, "Sum_"+tmp.m_name); //|A_i|^2
			newTree->createNode("NormFlatte_"+tmp.m_name, flattePhspStrat, "AbsVal_"+tmp.m_name, toyPhspSample.nEvents); //BW
			//Flatte (Normalization)
			newTree->createLeaf("m0_"+tmp.m_name, params.GetDoubleParameter("m0_"+name), "NormFlatte_"+tmp.m_name); //m0
			newTree->createLeaf("m23sq_phsp", m23sq_phsp, "NormFlatte_"+tmp.m_name); //ma
			newTree->createLeaf("m13sq_phsp", m13sq_phsp, "NormFlatte_"+tmp.m_name); //mb
			newTree->createLeaf("m12sq_phsp", m12sq_phsp, "NormFlatte_"+tmp.m_name); //mc
			newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "NormFlatte_"+tmp.m_name); //subSysFlag
			newTree->createLeaf("spin_"+tmp.m_name, tmp.m_spin, "NormFlatte_"+tmp.m_name); //spin
			newTree->createLeaf("d_"+tmp.m_name,  params.GetDoubleParameter("d_"+name), "NormFlatte_"+tmp.m_name); //d
			newTree->createLeaf("mHiddenA_"+tmp.m_name, \
					PhysConst::instance()->getMass(tmp.m_g2_part1), "NormFlatte_"+tmp.m_name);
			newTree->createLeaf("mHiddenB_"+tmp.m_name, \
					PhysConst::instance()->getMass(tmp.m_g2_part2), "NormFlatte_"+tmp.m_name);
//			newTree->createLeaf("g1_a_0", params.GetDoubleParameter("g1_a_0"), "NormFlatte_"+tmp.m_name);//use global parameter g1_a0 (asdfef)
			newTree->createLeaf("g1_a_0", params.GetDoubleParameter("g1_"+name), "NormFlatte_"+tmp.m_name);//use local parameter g1_a0
			newTree->createLeaf("g2_"+tmp.m_name, params.GetDoubleParameter("g2_"+name), "NormFlatte_"+tmp.m_name);
		} else {
			newTree->createLeaf("N_"+tmp.m_name, 1., "Flatte_"+tmp.m_name);
		}

		switch(subSys){
		case 3:{ //reso in sys of particles 1&2
			//newTree->createLeaf("mym_"+tmp.m_name, m12, "RelBW_"+tmp.m_name); //m
			newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "FlatteRes_"+tmp.m_name); //ma
			newTree->createLeaf("mb_"+tmp.m_name, kin->m2, "FlatteRes_"+tmp.m_name); //mb
			if(_normStyle!=normStyle::none){
				newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "NormFlatte_"+tmp.m_name); //ma
				newTree->createLeaf("mb_"+tmp.m_name, kin->m2, "NormFlatte_"+tmp.m_name); //mb
			}
			break;
		}
		case 4:{ //reso in sys of particles 1&3
			//newTree->createLeaf("mym_"+tmp.m_name, m13, "FlatteRes_"+tmp.m_name); //m
			newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "FlatteRes_"+tmp.m_name); //ma
			newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "FlatteRes_"+tmp.m_name); //mb
			if(_normStyle!=normStyle::none){
				newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "NormFlatte_"+tmp.m_name); //ma
				newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "NormFlatte_"+tmp.m_name); //mb
			}
			break;
		}
		case 5:{ //reso in sys of particles 2&3
			//newTree->createLeaf("mym_"+tmp.m_name, m23, "FlatteRes_"+tmp.m_name); //m
			newTree->createLeaf("ma_"+tmp.m_name, kin->m2, "FlatteRes_"+tmp.m_name); //ma
			newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "FlatteRes_"+tmp.m_name); //mb
			if(_normStyle!=normStyle::none){
				newTree->createLeaf("ma_"+tmp.m_name, kin->m2, "NormFlatte_"+tmp.m_name); //ma
				newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "NormFlatte_"+tmp.m_name); //mb
			}
			break;
		}
		default:{
			BOOST_LOG_TRIVIAL(error)<<"AmpSumIntensity::setupBasicTree(): Subsys not found!!";
		}
		}
	}
	BOOST_LOG_TRIVIAL(error)<<"AmpSumIntensity::setupBasicTree(): tree constructed!!";
	return newTree;
}

double AmpSumIntensity::getMaxVal(std::shared_ptr<Generator> gen){
	if(!_calcMaxFcnVal) calcMaxVal(gen);
	return _maxFcnVal;
}

double AmpSumIntensity::getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen){
	calcMaxVal(par,gen);
	return _maxFcnVal;
}

void AmpSumIntensity::calcMaxVal(ParameterList& par, std::shared_ptr<Generator> gen){
	setParameterList(par);
	return calcMaxVal(gen);
}

void AmpSumIntensity::calcMaxVal(std::shared_ptr<Generator> gen){
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());

	double maxM23=-999; double maxM13=-999; double maxVal=0;
	for(unsigned int i=0; i<_nCalls; i++){
		double m23sq=gen->getUniform()*(kin->m23_sq_max-kin->m23_sq_min)+kin->m23_sq_min;
		double m13sq=gen->getUniform()*(kin->m13_sq_max-kin->m13_sq_min)+kin->m13_sq_min;
		//		dataPoint point; point.setVal("m13sq",m13sq); point.setVal("m23sq",m23sq);
		dataPoint point; point.setVal(1,m13sq); point.setVal(0,m23sq);
		if( !kin->isWithinPhsp(point) ) { if(i>0) i--; continue; }//only integrate over phase space
		ParameterList res = intensity(point);
		double intens = *res.GetDoubleParameter(0);
		if(intens>maxVal) {
			maxM23=m23sq; maxM13=m13sq;
			maxVal=intens;
		}
	}
	_maxFcnVal=maxVal;
	_calcMaxFcnVal=1;
	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity::calcMaxVal() calculated maximum of amplitude: "
			<<_maxFcnVal<<" at m23sq="<<maxM23<<"/m13sq="<<maxM13;
	return ;
}
double AmpSumIntensity::normReso(std::shared_ptr<AmpAbsDynamicalFunction> amp){
	double norm;
	if(_normStyle==none) norm=amp->GetNormalization();
	else if(_normStyle==one) norm = sqrt(amp->integral(_nCalls));
	BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::normRes Normalization constant for "
			<<amp->GetName()<<": "<<1.0/norm;
	return norm;
}
/* Calculation amplitude integral (excluding efficiency) */
double AmpSumIntensity::evaluate(double x[], size_t dim) {
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//set data point: we assume that x[0]=m13 and x[1]=m23
	dataPoint point; point.setVal(1,x[0]); point.setVal(0,x[1]);
	//	double m12sq = kin->getThirdVariableSq(x[0],x[1]);
	if( !kin->isWithinPhsp(point) ) return 0;//only integrate over phase space
	ParameterList res = intensityNoEff(point);
	double intens = *res.GetDoubleParameter(0);
	return intens;
}
double evalWrapperAmpSumIntensity(double* x, size_t dim, void* param) {
	/* We need a wrapper here because intensity() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	return static_cast<AmpSumIntensity*>(param)->evaluate(x,dim);
};
const double AmpSumIntensity::integral(ParameterList& par){
	setParameterList(par);
	return integral();
}
const double AmpSumIntensity::integral(){
	/* Integration functionality was tested with a model with only one normalized amplitude.
	 * The integration result is equal to the amplitude coefficient^2.
	 */
	size_t dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	//	double xLimit_low[2] = {0,0};
	//	double xLimit_high[2] = {10,10};
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function G = {&evalWrapperAmpSumIntensity,dim, const_cast<AmpSumIntensity*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&G, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity::integrate() Integration result for amplitude sum: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}
/* Calculation amplitude normalization (including efficiency) */
double AmpSumIntensity::evaluateEff(double x[], size_t dim) {
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//set data point: we assume that x[0]=m13 and x[1]=m23
	dataPoint point; point.setVal(1,x[0]); point.setVal(0,x[1]);
	if( !kin->isWithinPhsp(point) ) return 0;//only integrate over phase space
	ParameterList res = intensity(point);
	double intens = *res.GetDoubleParameter(0);
	return intens;
}
double evalWrapperAmpSumIntensityEff(double* x, size_t dim, void* param) {
	/* We need a wrapper here because intensity() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	return static_cast<AmpSumIntensity*>(param)->evaluateEff(x,dim);
};

const double AmpSumIntensity::normalization(ParameterList& par){
	setParameterList(par);
	return integral();
}
const double AmpSumIntensity::normalization(){
	/* Integration functionality was tested with a model with only one normalized amplitude.
	 * The integration result is equal to the amplitude coefficient^2.
	 */
	size_t dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	//	double xLimit_low[2] = {0,0};
	//	double xLimit_high[2] = {10,10};
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function G = {&evalWrapperAmpSumIntensityEff,dim, const_cast<AmpSumIntensity*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&G, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity::normalization() Integration result for amplitude sum: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

std::complex<double> AmpSumIntensity::getFirstBW(dataPoint& point, ParameterList& par){
	setParameterList(par);
	return totAmp.getFirstBW(point);
}

std::complex<double> AmpSumIntensity::getFirstReso(dataPoint& point, ParameterList& par){
	setParameterList(par);
	return totAmp.getFirstReso(point);
}

std::complex<double> AmpSumIntensity::getFirstAmp(dataPoint& point, ParameterList& par){
	setParameterList(par);
	return totAmp.getFirstAmp(point);
}

const double AmpSumIntensity::sliceIntensity(dataPoint& dataP, ParameterList& par,std::complex<double>* reso, unsigned int nResos){
	setParameterList(par);
	//  dataPoint dataP; dataP.setVal("m23sq",point[0]); dataP.setVal("m13sq",point[1]);
	//dataPoint dataP; dataP.setVal(0,point[0]); dataP.setVal(1,point[1]);

	double AMPpdf=0;
	if(Kinematics::instance()->isWithinPhsp(dataP)) AMPpdf = totAmp.evaluateSlice(dataP, reso, nResos,5);
	if(AMPpdf!=AMPpdf){
		BOOST_LOG_TRIVIAL(error)<<"Error AmpSumIntensity: Intensity is not a number!!";
		AMPpdf = 0;
	}
	double eff=eff_->evaluate(dataP);
	return AMPpdf*eff;
}

const ParameterList& AmpSumIntensity::intensity(std::vector<double> point, ParameterList& par){
	setParameterList(par);
	dataPoint dataP; dataP.setVal(0,point[0]); dataP.setVal(1,point[1]);
	return intensity(dataP);
}
const ParameterList& AmpSumIntensity::intensity(dataPoint& point, ParameterList& par){
	setParameterList(par);
	return intensity(point);
}
const ParameterList& AmpSumIntensity::intensityNoEff(dataPoint& point){
	double AMPpdf=0;
	if(Kinematics::instance()->isWithinPhsp(point)) AMPpdf = totAmp.evaluate(point);

	if(AMPpdf!=AMPpdf){
		BOOST_LOG_TRIVIAL(error)<<"Error AmpSumIntensity: Intensity is not a number!!";
		AMPpdf = 0;
	}
	result.SetParameterValue(0,AMPpdf);
	return result;
}
const ParameterList& AmpSumIntensity::intensity(dataPoint& point){
	double AMPpdf=0;
	if(Kinematics::instance()->isWithinPhsp(point)) AMPpdf = totAmp.evaluate(point);

	if(AMPpdf!=AMPpdf){
		BOOST_LOG_TRIVIAL(error)<<"Error AmpSumIntensity: Intensity is not a number!!";
		AMPpdf = 0;
	}
	double eff=eff_->evaluate(point);
	result.SetParameterValue(0,AMPpdf*eff);
	return result;

}
void AmpSumIntensity::copyParameterList(ParameterList& outPar){
	if(outPar.GetNParameter())
		throw std::runtime_error("copyParameterList(): input list not empty");
	for(unsigned int i=0; i<params.GetNDouble();i++)
		outPar.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter(*(params.GetDoubleParameter(i)))) );

	return;
}

void AmpSumIntensity::setParameterList(ParameterList& par){
	//parameters varied by Minimization algorithm
	if(par.GetNDouble()!=params.GetNDouble())
		throw std::runtime_error("setParameterList(): size of parameter lists don't match");
	for(unsigned int i=0; i<params.GetNDouble(); i++){
		std::shared_ptr<DoubleParameter> p = params.GetDoubleParameter(i);
		if(!p->IsFixed()){
			p->SetValue(par.GetDoubleParameter(i)->GetValue());
			p->SetError(par.GetDoubleParameter(i)->GetError());
		}
	}
	return;
}
const bool AmpSumIntensity::fillStartParVec(ParameterList& outPar){
	outPar = ParameterList(params);
	return true;
}

void AmpSumIntensity::printAmps(){
	std::stringstream outStr;
	outStr<<"AmpSumIntensity: Printing amplitudes with current(!) set of parameters:\n";
	unsigned int n=0;
	for(unsigned int i=0; i<params.GetNDouble(); i++){
		std::shared_ptr<DoubleParameter> p = params.GetDoubleParameter(i);
		std::string tmp = p->GetName();
		std::size_t found = tmp.find("mag");
		if(found!=std::string::npos){
			outStr<<"-------- "<<namer[n]<<" ---------\n";
			n++;
		}
		outStr<<p->GetName()<<" = "<<p->GetValue();
		if(p->HasError())
			outStr<<"+-"<<p->GetError()->GetError();
		if(p->HasBounds())
			outStr<<" ["<<p->GetMinValue()<<";"<<p->GetMaxValue()<<"]";
		if(p->IsFixed())
			outStr<<" FIXED";
		outStr<<"\n";
	}

	BOOST_LOG_TRIVIAL(info)<<outStr.str();
	return;
}
void AmpSumIntensity::printFractions(){
	std::stringstream outStr;
	outStr<<"Fit fractions for all amplitudes: \n";
	double sumFrac=0;
	for(unsigned int i=0;i<totAmp.getNAmps();i++){
		double frac = getFraction(i);
		sumFrac+=frac;
		outStr<<std::setw(10)<<totAmp.getAmpName(i)<<":    "<<frac<<"\n";
		//		if(!(i==totAmp.getNAmps()-1)) outStr << "\n";
	}
	outStr<<std::setw(10)<<" "<<"    ==========\n";
	outStr<<std::setw(10)<<" "<<"     "<<sumFrac;
	BOOST_LOG_TRIVIAL(info)<<outStr.str();
	return;
}

double AmpSumIntensity::getIntValue(std::string var1, double min1, double max1, std::string var2, double min2, double max2){
	/*
	 * Integrates in \var1 from \min1 to \max1 and in \var2 from \min2 to \max2.
	 * Is intended to be used for calculation of bin content.
	 */
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double _min1 = min1;
	double _min2 = min2;
	double _max1 = max1;
	double _max2 = max2;
	//if(_min1==0) _min1 = kin->GetMin(var1);
	//if(_max1==0) _min1 = kin->GetMax(var1);
	if(_min2==0) _min2 = kin->getMin(var2);
	if(_max2==0) _max2 = kin->getMax(var2);
	unsigned int dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	double xLimit_low[2];
	double xLimit_high[2];

	if(var1 == "m13sq" && var2 == "m23sq"){
		xLimit_low[0] = _min1;
		xLimit_low[1] = _min2;
		xLimit_high[0] = _max1;
		xLimit_high[1] = _max2;
	}else if(var1 == "m23sq" && var2 == "m13sq"){
		xLimit_low[0] = _min2;
		xLimit_low[1] = _min1;
		xLimit_high[0] = _max2;
		xLimit_high[1] = _max1;
	} else {
		BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::getIntValue() wrong variables specified!";
		return -999;
	}

	size_t calls = 5000;
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function G = {&evalWrapperAmpSumIntensity,dim, const_cast<AmpSumIntensity*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&G, xLimit_low, xLimit_high, 2, calls, r,s,&res, &err);
	gsl_monte_vegas_free(s);

	return res;
}

/* OBSOLETE SECTION ONLY FOR TESTING */
//std::shared_ptr<FunctionTree> AmpSumIntensity::functionTree(allMasses& theMasses, allMasses& toyPhspSample) {
//	if(myTree) return myTree;
//	setupTree(theMasses, toyPhspSample, "data");
//	return myTree;
//}
//std::shared_ptr<FunctionTree> AmpSumIntensity::phspTree(allMasses& accPhspSample, allMasses& toyPhspSample) {
//	if(myPhspTree) return myPhspTree;
//	setupTree(accPhspSample,toyPhspSample,"normAcc");
//
//	return myPhspTree;
//}
//std::shared_ptr<FunctionTree> AmpSumIntensity::phspTree(allMasses& toyPhspSample) {
//	if(myPhspTree) return myPhspTree;
//	allMasses dummyMass;
//	setupTree(toyPhspSample,dummyMass,"norm");
//
//	return myPhspTree;
//}
//void AmpSumIntensity::setupTree(allMasses& theMasses, allMasses& toyPhspSample, std::string opt){
//	BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::setupTree() generating new tree!";
//	if(theMasses.nEvents==0){
//		BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity: sample empty!";
//		return;
//	}
//	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
//	_dpArea = kin->getPhspVolume();
//
//	if( opt == "data" ){
//		BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::setupTree() setting up data Tree";
//	} else if( opt == "norm" ){
//		if( toyPhspSample.nEvents!=0 ) throw(std::logic_error("Error in setupTree()! for the "
//				"normalization method the second sample should be size zero!"));
//		toyPhspSample = theMasses;
//		BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::setupTree() setting up normalization tree, "
//				"using efficiency corrected toy sample!";
//	} else if( opt == "normAcc" ){
//		BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::setupTree() setting up normalization tree, "
//				"using sample of accepted phsp events for efficiency correction!";
//	} else throw(std::logic_error("Error in setupTree()! Wrong option!"));
//
//	//------------Setup Tree---------------------
//	std::shared_ptr<FunctionTree> newTree = std::shared_ptr<FunctionTree>(new FunctionTree());
//	//------------Setup Tree Pars---------------------
//	std::shared_ptr<MultiDouble> m23sq = std::shared_ptr<MultiDouble>( new MultiDouble("m23sq",theMasses.masses_sq.at( std::make_pair(2,3) )) );
//	std::shared_ptr<MultiDouble> m13sq = std::shared_ptr<MultiDouble>( new MultiDouble("m13sq",theMasses.masses_sq.at( std::make_pair(1,3) )) );
//	std::shared_ptr<MultiDouble> m12sq = std::shared_ptr<MultiDouble>( new MultiDouble("m12sq",theMasses.masses_sq.at( std::make_pair(1,2) )) );
//	std::shared_ptr<MultiDouble> eff = std::shared_ptr<MultiDouble>( new MultiDouble("eff",theMasses.eff) ); //only needed for opt == "norm"
//	std::shared_ptr<MultiDouble> weight = std::shared_ptr<MultiDouble>( new MultiDouble("weight",theMasses.weight) );//only needed for opt == "data"
//	std::shared_ptr<MultiDouble> m23sq_phsp = std::shared_ptr<MultiDouble>( new MultiDouble("m23sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(2,3) )) );
//	std::shared_ptr<MultiDouble> m13sq_phsp = std::shared_ptr<MultiDouble>( new MultiDouble("m13sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(1,3) )) );
//	std::shared_ptr<MultiDouble> m12sq_phsp = std::shared_ptr<MultiDouble>( new MultiDouble("m12sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(1,2) )) );
//
//	//----Strategies needed
//	std::shared_ptr<MultAll> mmultStrat = std::shared_ptr<MultAll>(new MultAll(ParType::MCOMPLEX));
//	std::shared_ptr<MultAll> mmultDStrat = std::shared_ptr<MultAll>(new MultAll(ParType::MDOUBLE));
//	std::shared_ptr<AddAll> maddStrat = std::shared_ptr<AddAll>(new AddAll(ParType::MCOMPLEX));
//	std::shared_ptr<AbsSquare> msqStrat = std::shared_ptr<AbsSquare>(new AbsSquare(ParType::MDOUBLE));
//	std::shared_ptr<LogOf> mlogStrat = std::shared_ptr<LogOf>(new LogOf(ParType::MDOUBLE));
//	std::shared_ptr<MultAll> multStrat = std::shared_ptr<MultAll>(new MultAll(ParType::COMPLEX));
//	std::shared_ptr<MultAll> multDStrat = std::shared_ptr<MultAll>(new MultAll(ParType::DOUBLE));
//	std::shared_ptr<AddAll> addStrat = std::shared_ptr<AddAll>(new AddAll(ParType::DOUBLE));
//	std::shared_ptr<AddAll> addComplexStrat = std::shared_ptr<AddAll>(new AddAll(ParType::COMPLEX));
//	std::shared_ptr<AbsSquare> sqStrat = std::shared_ptr<AbsSquare>(new AbsSquare(ParType::DOUBLE));
//	std::shared_ptr<LogOf> logStrat = std::shared_ptr<LogOf>(new LogOf(ParType::DOUBLE));
//	std::shared_ptr<Complexify> complStrat = std::shared_ptr<Complexify>(new Complexify(ParType::COMPLEX));
//	std::shared_ptr<Inverse> invStrat = std::shared_ptr<Inverse>(new Inverse(ParType::DOUBLE));
//	std::shared_ptr<SquareRoot> sqRootStrat = std::shared_ptr<SquareRoot>(new SquareRoot(ParType::DOUBLE));
//
//	newTree->createHead("LH", addStrat); //Sum up all events, collapse multia
//
//	if( opt == "data" ){ //Data: EvtSum of log of Intens needed. Efficiency drops out in LH!
//		newTree->createNode("weightLog", mmultDStrat, "LH", theMasses.nEvents, false); //w_i * log( I_i )
//		newTree->createLeaf("weight", weight, "weightLog");
//		newTree->createNode("Log", mlogStrat, "weightLog", theMasses.nEvents, false); //log of amp, at each point
//		newTree->createNode("Intens", msqStrat, "Log", theMasses.nEvents, false); //I=A^2, at each point
//		//newTree->createNode("AmplitudeEff", mmultStrat, "Intens", theMasses.nEvents, false); //Sum of resonances * efficiency
//		//newTree->createLeaf("eff", eff, "AmplitudeEff"); //efficiency
//		//newTree->createNode("Amplitude", maddStrat, "AmplitudeEff", theMasses.nEvents, false); //Sum of resonances, at each point
//		newTree->createNode("Amplitude", maddStrat, "Intens", theMasses.nEvents, false); //Sum of resonances, at each point
//	} else if( opt == "norm" ){ //norm tree: sum of intensities, event efficiencies from toyPhspSample
//		newTree->createNode("Intens", msqStrat, "LH", theMasses.nEvents, false); //I=A^2, at each point
//		newTree->createNode("AmplitudeEff", mmultStrat, "Intens", theMasses.nEvents, false); //Sum of resonances, at each point
//		newTree->createLeaf("eff", eff, "AmplitudeEff"); //efficiency
//		newTree->createNode("Amplitude", maddStrat, "AmplitudeEff", theMasses.nEvents, false); //Sum of resonances, at each point
//	} else if( opt == "normAcc" ){ //norm tree: sum of intensities, theMasses should be a sample of accepted events
//		newTree->createNode("Intens", msqStrat, "LH", theMasses.nEvents, false); //I=A^2, at each point
//		newTree->createNode("Amplitude", maddStrat, "Intens", theMasses.nEvents, false); //Sum of resonances, at each point
//	} else throw(std::logic_error("Error in setupTree()! Wrong option!"));
//
////	std::vector<std::shared_ptr<DoubleParameter> >::iterator paramsPos = params.begin()+1;
////	std::vector<std::shared_ptr<DoubleParameter> >::iterator paramsPos = params.begin()+2;//use global parameter for g1
//
//	//----Add Resonances
//	for(std::vector<Resonance>::iterator reso=ampSetup.getResonances().begin(); reso!=ampSetup.getResonances().end(); reso++){
//		Resonance tmp = (*reso);
//		if(!tmp.m_enable) continue;
//		BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::setupTree() adding "<<tmp.m_name<<" (BW) to tree.";
//
//		//----Add Nodes
//		std::shared_ptr<BreitWignerStrategy> rbwStrat = std::shared_ptr<BreitWignerStrategy>(new BreitWignerStrategy(tmp.m_name,ParType::MCOMPLEX));
//		std::shared_ptr<BreitWignerPhspStrategy> rbwPhspStrat = std::shared_ptr<BreitWignerPhspStrategy>(new BreitWignerPhspStrategy(tmp.m_name,ParType::MCOMPLEX));
//		std::shared_ptr<WignerDStrategy> angdStrat = std::shared_ptr<WignerDStrategy>(new WignerDStrategy(tmp.m_name,ParType::MDOUBLE));
//		//		std::shared_ptr<WignerDphspStrategy> angdPhspStrat = std::shared_ptr<WignerDphspStrategy>(new WignerDphspStrategy(tmp.m_name,ParType::MDOUBLE));
//		unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;
//		newTree->createNode("Reso_"+tmp.m_name, mmultStrat, "Amplitude", theMasses.nEvents); //Reso=BW*C_*AD*N_
//		newTree->createNode("BW_"+tmp.m_name, mmultStrat , "Reso_"+tmp.m_name, theMasses.nEvents); //BW
//		newTree->createNode("RelBW_"+tmp.m_name, rbwStrat, "BW_"+tmp.m_name, theMasses.nEvents); //BW
//		newTree->createNode("C_"+tmp.m_name, complStrat, "Reso_"+tmp.m_name); //m0c
//		newTree->createLeaf("Intens_"+tmp.m_name, *paramsPos, "C_"+tmp.m_name); //r
//		newTree->createLeaf("Phase_"+tmp.m_name, *(paramsPos+1), "C_"+tmp.m_name); //phi
//		newTree->createNode("AngD_"+tmp.m_name, angdStrat, "Reso_"+tmp.m_name, theMasses.nEvents); //AD
//
//		//Breit-Wigner
//		newTree->createLeaf("m0_"+tmp.m_name, *(paramsPos+2), "RelBW_"+tmp.m_name); //m0
//		newTree->createLeaf("m23sq", m23sq, "RelBW_"+tmp.m_name); //ma
//		newTree->createLeaf("m13sq", m13sq, "RelBW_"+tmp.m_name); //mb
//		newTree->createLeaf("m12sq", m12sq, "RelBW_"+tmp.m_name); //mc
//		newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "RelBW_"+tmp.m_name); //subSysFlag
//		newTree->createLeaf("spin_"+tmp.m_name, tmp.m_spin, "RelBW_"+tmp.m_name); //spin
//		newTree->createLeaf("d_"+tmp.m_name, *(paramsPos+4) , "RelBW_"+tmp.m_name); //d
////		newTree->createLeaf("d_"+tmp.m_name, tmp.m_mesonRadius , "RelBW_"+tmp.m_name); //d
//		newTree->createLeaf("width_"+tmp.m_name, *(paramsPos+3), "RelBW_"+tmp.m_name); //resWidth
//		//Angular distribution
//		newTree->createLeaf("m23sq", m23sq, "AngD_"+tmp.m_name); //ma
//		newTree->createLeaf("m13sq", m13sq, "AngD_"+tmp.m_name); //mb
//		newTree->createLeaf("m12sq", m12sq, "AngD_"+tmp.m_name); //mc
//		newTree->createLeaf("M", kin->M, "AngD_"+tmp.m_name); //M
//		newTree->createLeaf("m1", kin->m1, "AngD_"+tmp.m_name); //m1
//		newTree->createLeaf("m2", kin->m2, "AngD_"+tmp.m_name); //m2
//		newTree->createLeaf("m3", kin->m3, "AngD_"+tmp.m_name); //m3
//		newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "AngD_"+tmp.m_name); //subSysFlag
//		newTree->createLeaf("spin_"+tmp.m_name,tmp.m_spin, "AngD_"+tmp.m_name); //spin
//		newTree->createLeaf("m_"+tmp.m_name, tmp.m_m, "AngD_"+tmp.m_name); //OutSpin 1
//		newTree->createLeaf("n_"+tmp.m_name, tmp.m_n, "AngD_"+tmp.m_name); //OutSpin 2
//
//		//adding nodes and leafs for calculation of normalization
//		if(_normStyle==normStyle::none){
//			newTree->createLeaf("N_"+tmp.m_name, 1., "BW_"+tmp.m_name);
//		}else{
//			//Normalization parameter for dynamical amplitude
//			newTree->createNode("N_"+tmp.m_name, sqRootStrat, "BW_"+tmp.m_name); //N = sqrt(NSq)
//			newTree->createNode("NSq_"+tmp.m_name, multDStrat, "N_"+tmp.m_name); //NSq = N_phspMC * 1/PhspVolume * 1/Sum(|A|^2)
//			newTree->createLeaf("PhspSize_"+tmp.m_name, toyPhspSample.nEvents, "NSq_"+tmp.m_name); // N_phspMC
//			newTree->createLeaf("PhspVolume_"+tmp.m_name, 1/_dpArea, "NSq_"+tmp.m_name); // 1/PhspVolume
//			newTree->createNode("InvSum_"+tmp.m_name, invStrat, "NSq_"+tmp.m_name); //1/Sum(|A|^2)
//			newTree->createNode("Sum_"+tmp.m_name, addStrat, "InvSum_"+tmp.m_name); //Sum(|A|^2)
//			newTree->createNode("AbsVal_"+tmp.m_name, msqStrat, "Sum_"+tmp.m_name); //|A_i|^2
//			//Breit-Wigner (Normalization)
//			newTree->createNode("NormBW_"+tmp.m_name, rbwPhspStrat, "AbsVal_"+tmp.m_name, toyPhspSample.nEvents); //BW
//			newTree->createLeaf("m0_"+tmp.m_name, *(paramsPos+2), "NormBW_"+tmp.m_name); //m0
//			newTree->createLeaf("m23sq_phsp", m23sq_phsp, "NormBW_"+tmp.m_name); //ma
//			newTree->createLeaf("m13sq_phsp", m13sq_phsp, "NormBW_"+tmp.m_name); //mb
//			newTree->createLeaf("m12sq_phsp", m12sq_phsp, "NormBW_"+tmp.m_name); //mc
//			newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "NormBW_"+tmp.m_name); //subSysFlag
//			newTree->createLeaf("spin_"+tmp.m_name, tmp.m_spin, "NormBW_"+tmp.m_name); //spin
//			newTree->createLeaf("d_"+tmp.m_name, *(paramsPos+4) , "NormBW_"+tmp.m_name); //d
//			newTree->createLeaf("width_"+tmp.m_name, *(paramsPos+3), "NormBW_"+tmp.m_name); //resWidth
//		}
//		switch(subSys){
//		case 3:{ //reso in sys of particles 1&2
//			//newTree->createLeaf("mym_"+tmp.m_name, m12, "RelBW_"+tmp.m_name); //m
//			newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "RelBW_"+tmp.m_name); //ma
//			newTree->createLeaf("mb_"+tmp.m_name, kin->m2, "RelBW_"+tmp.m_name); //mb
//			if(_normStyle!=normStyle::none){
//				newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "NormBW_"+tmp.m_name); //ma
//				newTree->createLeaf("mb_"+tmp.m_name, kin->m2, "NormBW_"+tmp.m_name); //mb
//			}
//			break;
//		}
//		case 4:{ //reso in sys of particles 1&3
//			//newTree->createLeaf("mym_"+tmp.m_name, m13, "RelBW_"+tmp.m_name); //m
//			newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "RelBW_"+tmp.m_name); //ma
//			newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "RelBW_"+tmp.m_name); //mb
//			if(_normStyle!=normStyle::none){
//				newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "NormBW_"+tmp.m_name); //ma
//				newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "NormBW_"+tmp.m_name); //mb
//			}
//			break;
//		}
//		case 5:{ //reso in sys of particles 2&3
//			//newTree->createLeaf("mym_"+tmp.m_name, m23, "RelBW_"+tmp.m_name); //m
//			newTree->createLeaf("ma_"+tmp.m_name, kin->m2, "RelBW_"+tmp.m_name); //ma
//			newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "RelBW_"+tmp.m_name); //mb
//			if(_normStyle!=normStyle::none){
//				newTree->createLeaf("ma_"+tmp.m_name, kin->m2, "NormBW_"+tmp.m_name); //ma
//				newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "NormBW_"+tmp.m_name); //mb
//			}
//			break;
//		}
//		default:{
//			BOOST_LOG_TRIVIAL(error)<<"AmpSumIntensity::setupTree(): Subsys not found!!";
//		}
//		}
//		paramsPos += 5;
//	}// end loop over resonances
//
//	for(std::vector<ResonanceFlatte>::iterator reso=ampSetup.getResonancesFlatte().begin(); reso!=ampSetup.getResonancesFlatte().end(); reso++){
//		ResonanceFlatte tmp = (*reso);
//		if(!tmp.m_enable) continue;
//		BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::setupTree() adding "<<tmp.m_name<<" (Flatte) to tree.";
//
//		//----Add Nodes
//		std::shared_ptr<FlatteStrategy> flatteStrat = std::shared_ptr<FlatteStrategy>(new FlatteStrategy(tmp.m_name,ParType::MCOMPLEX));
//		std::shared_ptr<FlattePhspStrategy> flattePhspStrat = std::shared_ptr<FlattePhspStrategy>(new FlattePhspStrategy(tmp.m_name,ParType::MCOMPLEX));
//		std::shared_ptr<WignerDStrategy> angdStrat = std::shared_ptr<WignerDStrategy>(new WignerDStrategy(tmp.m_name,ParType::MDOUBLE));
//		//		std::shared_ptr<WignerDphspStrategy> angdPhspStrat = std::shared_ptr<WignerDphspStrategy>(new WignerDphspStrategy(tmp.m_name,ParType::MDOUBLE));
//		unsigned int subSys = tmp.m_daugtherA + tmp.m_daugtherB;
//		newTree->createNode("Reso_"+tmp.m_name, mmultStrat, "Amplitude", theMasses.nEvents); //Reso=BW*C_*AD*N_
//		newTree->createNode("Flatte_"+tmp.m_name, mmultStrat , "Reso_"+tmp.m_name, theMasses.nEvents); //BW
//		newTree->createNode("FlatteRes_"+tmp.m_name, flatteStrat, "Flatte_"+tmp.m_name, theMasses.nEvents); //BW
//		newTree->createNode("C_"+tmp.m_name, complStrat, "Reso_"+tmp.m_name); //c
//		newTree->createLeaf("Intens_"+tmp.m_name, *paramsPos, "C_"+tmp.m_name); //r
//		newTree->createLeaf("Phase_"+tmp.m_name, *(paramsPos+1), "C_"+tmp.m_name); //phi
//		newTree->createNode("AngD_"+tmp.m_name, angdStrat, "Reso_"+tmp.m_name, theMasses.nEvents); //AD
//
//		//Flatte
//		newTree->createLeaf("m0_"+tmp.m_name, *(paramsPos+2), "FlatteRes_"+tmp.m_name); //m0
//		newTree->createLeaf("m23sq", m23sq, "FlatteRes_"+tmp.m_name); //ma
//		newTree->createLeaf("m13sq", m13sq, "FlatteRes_"+tmp.m_name); //mb
//		newTree->createLeaf("m12sq", m12sq, "FlatteRes_"+tmp.m_name); //mc
//		newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "FlatteRes_"+tmp.m_name); //subSysFlag
//		newTree->createLeaf("spin_"+tmp.m_name, tmp.m_spin, "FlatteRes_"+tmp.m_name); //spin
//		newTree->createLeaf("d_"+tmp.m_name, *(paramsPos+5) , "FlatteRes_"+tmp.m_name); //d
//		newTree->createLeaf("mHiddenA_"+tmp.m_name, \
//				PhysConst::instance()->getMass(tmp.m_g2_part1), "FlatteRes_"+tmp.m_name);
//		newTree->createLeaf("mHiddenB_"+tmp.m_name, \
//				PhysConst::instance()->getMass(tmp.m_g2_part2), "FlatteRes_"+tmp.m_name);
////		newTree->createLeaf("g1_"+tmp.m_name, *(paramsPos+3), "FlatteRes_"+tmp.m_name);
//		newTree->createLeaf("g1_"+tmp.m_name, params.at(1), "FlatteRes_"+tmp.m_name);
//		newTree->createLeaf("g2_"+tmp.m_name, *(paramsPos+4), "FlatteRes_"+tmp.m_name);
//		//Angular distribution
//		newTree->createLeaf("m23sq", m23sq, "AngD_"+tmp.m_name); //ma
//		newTree->createLeaf("m13sq", m13sq, "AngD_"+tmp.m_name); //mb
//		newTree->createLeaf("m12sq", m12sq, "AngD_"+tmp.m_name); //mc
//		newTree->createLeaf("M", kin->M, "AngD_"+tmp.m_name); //M
//		newTree->createLeaf("m1", kin->m1, "AngD_"+tmp.m_name); //m1
//		newTree->createLeaf("m2", kin->m2, "AngD_"+tmp.m_name); //m2
//		newTree->createLeaf("m3", kin->m3, "AngD_"+tmp.m_name); //m3
//		newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "AngD_"+tmp.m_name); //subSysFlag
//		newTree->createLeaf("spin_"+tmp.m_name,tmp.m_spin, "AngD_"+tmp.m_name); //spin
//		newTree->createLeaf("m_"+tmp.m_name, tmp.m_m, "AngD_"+tmp.m_name); //OutSpin 1
//		newTree->createLeaf("n_"+tmp.m_name, tmp.m_n, "AngD_"+tmp.m_name); //OutSpin 2
//
//		//Normalization
//		if(_normStyle!=normStyle::none){
//			newTree->createNode("N_"+tmp.m_name, sqRootStrat, "Flatte_"+tmp.m_name); //N = sqrt(NSq)
//			newTree->createNode("NSq_"+tmp.m_name, multDStrat, "N_"+tmp.m_name); //NSq = N_phspMC * 1/PhspVolume * 1/Sum(|A|^2)
//			newTree->createLeaf("PhspSize_"+tmp.m_name, toyPhspSample.nEvents, "NSq_"+tmp.m_name); // N_phspMC
//			newTree->createLeaf("PhspVolume_"+tmp.m_name, 1/_dpArea, "NSq_"+tmp.m_name); // 1/PhspVolume
//			newTree->createNode("InvSum_"+tmp.m_name, invStrat, "NSq_"+tmp.m_name); //1/Sum(|A|^2)
//			newTree->createNode("Sum_"+tmp.m_name, addStrat, "InvSum_"+tmp.m_name); //Sum(|A|^2)
//			newTree->createNode("AbsVal_"+tmp.m_name, msqStrat, "Sum_"+tmp.m_name); //|A_i|^2
//			newTree->createNode("NormFlatte_"+tmp.m_name, flattePhspStrat, "AbsVal_"+tmp.m_name, toyPhspSample.nEvents); //BW
//			//Flatte (Normalization)
//			newTree->createLeaf("m0_"+tmp.m_name, *(paramsPos+2), "NormFlatte_"+tmp.m_name); //m0
//			newTree->createLeaf("m23sq_phsp", m23sq_phsp, "NormFlatte_"+tmp.m_name); //ma
//			newTree->createLeaf("m13sq_phsp", m13sq_phsp, "NormFlatte_"+tmp.m_name); //mb
//			newTree->createLeaf("m12sq_phsp", m12sq_phsp, "NormFlatte_"+tmp.m_name); //mc
//			newTree->createLeaf("subSysFlag_"+tmp.m_name, subSys, "NormFlatte_"+tmp.m_name); //subSysFlag
//			newTree->createLeaf("spin_"+tmp.m_name, tmp.m_spin, "NormFlatte_"+tmp.m_name); //spin
//			newTree->createLeaf("d_"+tmp.m_name,  *(paramsPos+5), "NormFlatte_"+tmp.m_name); //d
//			newTree->createLeaf("mHiddenA_"+tmp.m_name, \
//					PhysConst::instance()->getMass(tmp.m_g2_part1), "NormFlatte_"+tmp.m_name);
//			newTree->createLeaf("mHiddenB_"+tmp.m_name, \
//					PhysConst::instance()->getMass(tmp.m_g2_part2), "NormFlatte_"+tmp.m_name);
////			newTree->createLeaf("g1_"+tmp.m_name, *(paramsPos+3), "NormFlatte_"+tmp.m_name);
//			newTree->createLeaf("g1_"+tmp.m_name, params.at(1), "NormFlatte_"+tmp.m_name);
//			newTree->createLeaf("g2_"+tmp.m_name, *(paramsPos+4), "NormFlatte_"+tmp.m_name);
//		} else {
//			newTree->createLeaf("N_"+tmp.m_name, 1., "Flatte_"+tmp.m_name);
//		}
//
//		switch(subSys){
//		case 3:{ //reso in sys of particles 1&2
//			//newTree->createLeaf("mym_"+tmp.m_name, m12, "RelBW_"+tmp.m_name); //m
//			newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "FlatteRes_"+tmp.m_name); //ma
//			newTree->createLeaf("mb_"+tmp.m_name, kin->m2, "FlatteRes_"+tmp.m_name); //mb
//			if(_normStyle!=normStyle::none){
//				newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "NormFlatte_"+tmp.m_name); //ma
//				newTree->createLeaf("mb_"+tmp.m_name, kin->m2, "NormFlatte_"+tmp.m_name); //mb
//			}
//			break;
//		}
//		case 4:{ //reso in sys of particles 1&3
//			//newTree->createLeaf("mym_"+tmp.m_name, m13, "FlatteRes_"+tmp.m_name); //m
//			newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "FlatteRes_"+tmp.m_name); //ma
//			newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "FlatteRes_"+tmp.m_name); //mb
//			if(_normStyle!=normStyle::none){
//				newTree->createLeaf("ma_"+tmp.m_name, kin->m1, "NormFlatte_"+tmp.m_name); //ma
//				newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "NormFlatte_"+tmp.m_name); //mb
//			}
//			break;
//		}
//		case 5:{ //reso in sys of particles 2&3
//			//newTree->createLeaf("mym_"+tmp.m_name, m23, "FlatteRes_"+tmp.m_name); //m
//			newTree->createLeaf("ma_"+tmp.m_name, kin->m2, "FlatteRes_"+tmp.m_name); //ma
//			newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "FlatteRes_"+tmp.m_name); //mb
//			if(_normStyle!=normStyle::none){
//				newTree->createLeaf("ma_"+tmp.m_name, kin->m2, "NormFlatte_"+tmp.m_name); //ma
//				newTree->createLeaf("mb_"+tmp.m_name, kin->m3, "NormFlatte_"+tmp.m_name); //mb
//			}
//			break;
//		}
//		default:{
//			BOOST_LOG_TRIVIAL(error)<<"AmpSumIntensity::setupTree(): Subsys not found!!";
//		}
//		}
//		paramsPos += 6;
//	}
//	if( opt == "data") myTree=newTree;
//	else if( opt == "norm" || opt == "normAcc") myPhspTree=newTree;
//	else throw(std::logic_error("Error in setupTree()! Wrong option!"));
//}
/* OBSOLETE SECTION ONLY FOR TESTING */
