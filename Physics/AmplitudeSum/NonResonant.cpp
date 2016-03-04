/*
 * NonResonant.cpp
 *
 *  Created on: Jan 13, 2015
 *      Author: weidenka
 */

#include "Physics/AmplitudeSum/NonResonant.hpp"


NonResonant::NonResonant(const char *name,
		std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
		int nCalls, normStyle nS) :
		AmpAbsDynamicalFunction(name,mag, phase,
				std::make_shared<DoubleParameter>("mass", 0.0), 1, 2, Spin(0), Spin(0), Spin(0),
				formFactorType::noFormFactor, nCalls, nS)
{
}

std::complex<double> NonResonant::dynamicalFunction(){
	return std::complex<double>(1,0);
}

//! Configure resonance from ptree
void NonResonant::Configure(boost::property_tree::ptree::value_type const& v,
		ParameterList& list)
{
	if( v.first != "NonResonant" ) throw BadConfig("");

	boost::property_tree::ptree pt = v.second;
	AmpAbsDynamicalFunction::Configure(v,list);

	return;
}

void NonResonant::Save(boost::property_tree::ptree &pt)
{
	boost::property_tree::ptree amp;
	AmpAbsDynamicalFunction::put(amp);
	pt.add_child("NonResonant", amp);
	return;
}

std::shared_ptr<FunctionTree> NonResonant::SetupTree(
		allMasses& theMasses,allMasses& toyPhspSample,std::string suffix){

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double phspVol = kin->getPhspVolume();
	BOOST_LOG_TRIVIAL(info) << "NonResonant::setupBasicTree() | "<<_name;
	//------------Setup Tree---------------------
	std::shared_ptr<FunctionTree> newTree(new FunctionTree());

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

	newTree->createHead("Reso_"+_name, mmultStrat, theMasses.nEvents); //Reso=C*N*nonReso
	newTree->createNode("C_"+_name, complStrat, "Reso_"+_name); //c=r*exp(phi)
	newTree->createLeaf("Intens_"+_name, _mag, "C_"+_name); //r
	newTree->createLeaf("Phase_"+_name, _phase, "C_"+_name); //phi

	std::shared_ptr<MultiComplex> unitVec(
			new MultiComplex("unit",std::vector<std::complex<double> >(
					theMasses.nEvents, std::complex<double>(1,0))) );

	newTree->createLeaf("NonRes_"+_name, unitVec, "Reso_"+_name); //nonReso
	//adding nodes and leafs for calculation of normalization
	if(_normStyle==normStyle::none){
		newTree->createLeaf("N_"+_name, 1., "Reso_"+_name);
	}else{
		newTree->createNode("N_"+_name, sqRootStrat, "Reso_"+_name); //N = sqrt(NSq)
		newTree->createNode("NSq_"+_name, multDStrat, "N_"+_name); //NSq = (PhspVolume/N_phspMC * Sum(|A|^2))^-1
		newTree->createLeaf("PhspSize_"+_name, toyPhspSample.nEvents, "NSq_"+_name); // N_phspMC
		newTree->createLeaf("PhspVolume_"+_name, 1/phspVol, "NSq_"+_name); // 1/PhspVolume
		newTree->createNode("InvSum_"+_name, invStrat, "NSq_"+_name); //1/Sum(|A|^2)
		newTree->createNode("Sum_"+_name, addStrat, "InvSum_"+_name); //Sum(|A|^2)
		newTree->createNode("AbsVal_"+_name, msqStrat, "Sum_"+_name); //|A_i|^2
		std::shared_ptr<MultiComplex> unitVec2(
				new MultiComplex("unit",std::vector<std::complex<double> >(
						toyPhspSample.nEvents, std::complex<double>(1,0))) );
		newTree->createLeaf("NormNonRes_"+_name, unitVec2, "AbsVal_"+_name); //BW
	}
	return newTree;
}
