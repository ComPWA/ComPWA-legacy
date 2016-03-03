//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//		Peter Weidenkaff - correct nominator, using dataPoint for data handling
//-------------------------------------------------------------------------------
//****************************************************************************
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Class for defining the relativistic Breit-Wigner resonance model, which
// includes the use of Blatt-Weisskopf barrier factors.

#include <cmath>
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include <stdlib.h>

AmpRelBreitWignerRes::AmpRelBreitWignerRes(const char *name,
		std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
		std::shared_ptr<DoubleParameter> mass, int part1, int part2,
		Spin spin, Spin m, Spin n,
		std::shared_ptr<DoubleParameter> width,
		std::shared_ptr<DoubleParameter> mesonRadius,
		std::shared_ptr<DoubleParameter> motherRadius,
		formFactorType type,
		int nCalls, normStyle nS) :
		AmpAbsDynamicalFunction(name, mag, phase, mass, part1, part2, spin, m, n,
				mesonRadius, motherRadius, type, nCalls, nS),
				_width(width)
{
	if( _width->GetValue() != tmp_width) {
		SetModified();
		tmp_width = _width->GetValue();
	}
	//setting default normalization
	GetNormalization();
}

AmpRelBreitWignerRes::~AmpRelBreitWignerRes() 
{
}

void AmpRelBreitWignerRes::Configure(
		boost::property_tree::ptree::value_type const& v,
		ParameterList& list)
{
	if( v.first != "BreitWigner" ) throw BadConfig("");

	boost::property_tree::ptree pt = v.second;
	AmpAbsDynamicalFunction::Configure(v,list);

	//Width (mandatory)
	auto tmp_width_fix = pt.get<bool>("width_fix",1);
	auto tmp_width_min = pt.get<double>("width_min",0.0);
	auto tmp_width_max = pt.get<double>("width_max",1.0);
	auto tmp_width_name = pt.get_optional<std::string>("width_name");
	if(!tmp_width_name){
		auto tmp_width= pt.get_optional<double>("width");
		if(!tmp_width)
			throw BadParameter("AmpRelBreitWignerRes::Configure() | "
					"width for "+_name+" not specified!");
		_width = std::shared_ptr<DoubleParameter>(
				new DoubleParameter(
						"width_"+_name,tmp_width.get(),
						tmp_width_min, tmp_width_max
				)
		);
		_width->FixParameter(tmp_width_fix);
		if(_enable) list.AddParameter(_width);
		_width_writeByName = 0;
	} else {
		try{
			_width = list.GetDoubleParameter(tmp_width_name.get());
			_width_writeByName = 1;
		} catch (BadParameter& ex){
			BOOST_LOG_TRIVIAL(error) <<"AmpRelBreitWignerRes::Configure() | "
					"Requesting parameter "<<tmp_width_name.get()<<" but"
							" was not found in parameter list. "
							"Quit since parameter is mandatory!";
			throw;
		}
	}

	initialize();
}

void AmpRelBreitWignerRes::Save(boost::property_tree::ptree &pt)
{
	boost::property_tree::ptree amp;
	AmpAbsDynamicalFunction::put(amp);
	if(_width_writeByName){
		amp.put("width_name", _width->GetName());
	} else {
		amp.put("width", _width->GetValue());
		amp.put("width_fix", _width->IsFixed());
		amp.put("width_min", _width->GetMinValue());
		amp.put("width_max", _width->GetMaxValue());
	}

	pt.add_child("BreitWigner", amp);
	return;
}

std::string AmpRelBreitWignerRes::to_str() const{
	std::string dynAmp = AmpAbsDynamicalFunction::to_str();
	std::stringstream str;
	str<<_width->to_str()<<std::endl;

	return dynAmp+str.str();
}

void AmpRelBreitWignerRes::CheckModified() {
	AmpAbsDynamicalFunction::CheckModified();
	if(_width->GetValue() != tmp_width){
		SetModified();
		tmp_width = _width->GetValue();
	}
	return;
}

std::complex<double> AmpRelBreitWignerRes::EvaluateAmp(dataPoint& point) {
	CheckModified(); //recalculate normalization ?

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double mSq = -999;
	switch(_subSys){
	case 3: mSq=kin->getThirdVariableSq(point.getVal(0),point.getVal(1)); break;
	case 4: mSq=point.getVal(1); break;
	case 5: mSq=point.getVal(0); break;
	}

	return dynamicalFunction(mSq,_mass->GetValue(),_ma,_mb,_width->GetValue(),_spin,
			_mesonRadius->GetValue(), _ffType);
}
std::complex<double> AmpRelBreitWignerRes::dynamicalFunction(double mSq, double mR,
		double ma, double mb, double width, unsigned int J, double mesonRadius,
		formFactorType ffType)
{
	std::complex<double> i(0,1);
	double sqrtS = sqrt(mSq);

	double barrier = Kinematics::FormFactor(sqrtS,ma,mb,J,mesonRadius, ffType)/
			Kinematics::FormFactor(mR,ma,mb,J,mesonRadius, ffType);

	std::complex<double> qTerm = std::pow(
			(Kinematics::phspFactor(sqrtS,ma,mb) / Kinematics::phspFactor(mR,ma,mb)) *mR/sqrtS,
			(2*J+ 1)
			);

	//Calculate coupling constant to final state
	std::complex<double> g_final = widthToCoupling(
			mSq,mR,width,ma,mb,J,mesonRadius,ffType
			);

	/*Coupling constant from production reaction. In case of a particle decay
	 * the production coupling doesn't depend in energy since the CM energy
	 * is in the (RC) system fixed to the mass of the decaying particle */
	double g_production = 1;

	std::complex<double> denom =
			std::complex<double>( mR*mR - mSq,0)
			+ (-1.0)*i*sqrtS*(width*qTerm*barrier);

	std::complex<double> result = g_final*g_production / denom;

	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"AmpRelBreitWignerRes::dynamicalFunction() | "
				<< barrier<<" "<<mR<<" "<<mSq <<" "<<ma<<" "<<mb<<std::endl;
		return 0;
	}
	return result;
}

std::shared_ptr<FunctionTree> AmpRelBreitWignerRes::SetupTree(
		allMasses& theMasses,allMasses& toyPhspSample,std::string suffix)
{
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double phspVol = kin->getPhspVolume();
	BOOST_LOG_TRIVIAL(info) << "AmpRelBreitWignerRes::setupBasicTree() | "
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
	std::shared_ptr<BreitWignerStrategy> rbwStrat(
			new BreitWignerStrategy(_name,ParType::MCOMPLEX) );
	std::shared_ptr<WignerDStrategy> angdStrat(
			new WignerDStrategy(_name,ParType::MDOUBLE) );

	//Reso=BW*C*AD*N
	newTree->createHead("Reso_"+_name, mmultStrat, theMasses.nEvents);
	newTree->createNode("C_"+_name, complStrat, "Reso_"+_name); //m0c
	newTree->createLeaf("Intens_"+_name, _mag, "C_"+_name); //r
	newTree->createLeaf("Phase_"+_name, _phase, "C_"+_name); //phi
	//Angular distribution
	newTree->insertTree(_wignerD.SetupTree(theMasses,suffix), "Reso_"+_name);

	//Breit-Wigner
	newTree->createNode("RelBW_"+_name, rbwStrat, "Reso_"+_name, theMasses.nEvents);
	newTree->createLeaf("mass", _mass, "RelBW_"+_name); //m0
	newTree->createLeaf("width", _width, "RelBW_"+_name); //resWidth
	newTree->createLeaf("spin", _spin, "RelBW_"+_name); //spin
	newTree->createLeaf("mesonRadius", _mesonRadius, "RelBW_"+_name); //d
	newTree->createLeaf("formFactorType", _ffType , "RelBW_"+_name); //d
	newTree->createLeaf("subSysFlag", _subSys, "RelBW_"+_name); //subSysFlag
	switch(_subSys){
	case 3:{ //reso in sys of particles 1&2
		newTree->createLeaf("ma", kin->m1, "RelBW_"+_name); //ma
		newTree->createLeaf("mb", kin->m2, "RelBW_"+_name); //mb
		break;
	}
	case 4:{ //reso in sys of particles 1&3
		newTree->createLeaf("ma", kin->m1, "RelBW_"+_name); //ma
		newTree->createLeaf("mb", kin->m3, "RelBW_"+_name); //mb
		break;
	}
	case 5:{ //reso in sys of particles 2&3
		newTree->createLeaf("ma", kin->m2, "RelBW_"+_name); //ma
		newTree->createLeaf("mb", kin->m3, "RelBW_"+_name); //mb
		break;
	}
	default:{
		BOOST_LOG_TRIVIAL(error)<<"AmpRelBreitWignerRes::setupBasicTree() | "
				"Subsys not found!";
	}
	}
	newTree->createLeaf("m12sq", m12sq, "RelBW_"+_name); //mc
	newTree->createLeaf("m13sq", m13sq, "RelBW_"+_name); //mb
	newTree->createLeaf("m23sq", m23sq, "RelBW_"+_name); //ma

	//Normalization
	if(_normStyle==normStyle::none) {
		newTree->createLeaf("N_"+_name, 1., "Reso_"+_name);
	} else {
		//Normalization parameter for dynamical amplitude
		newTree->createNode("N_"+_name, sqRootStrat, "Reso_"+_name); //N = sqrt(NSq)
		newTree->createNode("NSq_"+_name, multDStrat, "N_"+_name); //NSq = (PhspVolume/N_phspMC * Sum(|A|^2))^-1
		newTree->createLeaf("PhspSize_"+_name, toyPhspSample.nEvents, "NSq_"+_name); // N_phspMC
		newTree->createLeaf("PhspVolume_"+_name, 1/phspVol, "NSq_"+_name); // 1/PhspVolume
		newTree->createNode("InvSum_"+_name, invStrat, "NSq_"+_name); //1/Sum(|A|^2)
		newTree->createNode("Sum_"+_name, addStrat, "InvSum_"+_name); //Sum(|A|^2)
		newTree->createNode("AbsVal_"+_name, msqStrat, "Sum_"+_name); //|A_i|^2

		newTree->createNode("NormReso_"+_name, mmultStrat, "AbsVal_"+_name,
				toyPhspSample.nEvents);
		//Angular distribution (Normalization)
		newTree->insertTree(_wignerD.SetupTree(toyPhspSample,suffix),
				"NormReso_"+_name);
		//Breit-Wigner (Normalization)
		newTree->createNode("NormBW_"+_name, rbwStrat, "NormReso_"+_name,
				toyPhspSample.nEvents); //BW
		newTree->createLeaf("mass", _mass, "NormBW_"+_name); //m0
		newTree->createLeaf("width", _width, "NormBW_"+_name); //resWidth
		newTree->createLeaf("spin", _spin, "NormBW_"+_name); //spin
		newTree->createLeaf("mesonRadius", _mesonRadius, "NormBW_"+_name); //d
		newTree->createLeaf("formFactorType", _ffType , "NormBW_"+_name); //d
		newTree->createLeaf("subSysFlag", _subSys, "NormBW_"+_name); //subSysFlag
		switch(_subSys){
		case 3:{ //reso in sys of particles 1&2
			newTree->createLeaf("ma", kin->m1, "NormBW_"+_name); //ma
			newTree->createLeaf("mb", kin->m2, "NormBW_"+_name); //mb
			break;
		} case 4: { //reso in sys of particles 1&3
			newTree->createLeaf("ma", kin->m1, "NormBW_"+_name); //ma
			newTree->createLeaf("mb", kin->m3, "NormBW_"+_name); //mb
			break;
		} case 5: { //reso in sys of particles 2&3
			newTree->createLeaf("ma", kin->m2, "NormBW_"+_name); //ma
			newTree->createLeaf("mb", kin->m3, "NormBW_"+_name); //mb
			break;
		}
		default:{
			BOOST_LOG_TRIVIAL(error)<<"AmpSumIntensity::setupBasicTree(): "
					"Subsys not found!!";
		}
		}
		newTree->createLeaf("m12sq_phsp", m12sq_phsp, "NormBW_"+_name); //mc
		newTree->createLeaf("m13sq_phsp", m13sq_phsp, "NormBW_"+_name); //mb
		newTree->createLeaf("m23sq_phsp", m23sq_phsp, "NormBW_"+_name); //ma
	}

	BOOST_LOG_TRIVIAL(debug) << "AmpRelBreitWignerRes::setupBasicTree() | "
			"finished!";

	return newTree;
}

bool BreitWignerStrategy::execute(ParameterList& paras,
		std::shared_ptr<AbsParameter>& out)
{
	//Debug
	//	BOOST_LOG_TRIVIAL(debug) <<"BreitWignerStrategy::execute() | start";
	//	for( int i=0; i<paras.GetDoubleParameters().size(); ++i){
	//		BOOST_LOG_TRIVIAL(debug)<<paras.GetDoubleParameter(i)->GetName();
	//	}

	//Check parameter type
	if( checkType != out->type() )
		throw(WrongParType(std::string("Output Type ")
	+ ParNames[out->type()] + std::string(" conflicts expected type ")
	+ ParNames[checkType] + std::string(" of ")+name+" BW strat")
		);

	//Check size of parameter list
	if( paras.GetNDouble() != 8 && paras.GetNDouble() != 11)
		throw( BadParameter("BreitWignerStrategy::execute() | "
				"number of DoubleParameters does not match!")
		);

	double Gamma0, m0, d, ma, mb;
	unsigned int spin, subSys;
	int ffType;
	/** Get parameters from ParameterList:
	 * We use the same order of the parameters as was used during tree
	 * construction
	 */
	m0 = paras.GetDoubleParameter(0)->GetValue();
	Gamma0 = paras.GetDoubleParameter(1)->GetValue();
	spin = (unsigned int) paras.GetDoubleParameter(2)->GetValue();
	d = paras.GetDoubleParameter(3)->GetValue();
	ffType = paras.GetDoubleParameter(4)->GetValue();
	subSys = paras.GetDoubleParameter(5)->GetValue();
	ma = paras.GetDoubleParameter(6)->GetValue();
	mb = paras.GetDoubleParameter(7)->GetValue();

//	BOOST_LOG_TRIVIAL(debug) << "BreitWignerStrategy::execute() | mR="<<m0
//			<<" Gamma="<<Gamma0<<" spin="<<spin<<" radius="<<d<<" ffType="<<ffType
//			<<" subSys="<<subSys<<" ma="<<ma<<" mb="<<mb;

	//MultiDim output, must have multidim Paras in input
	if(checkType == ParType::MCOMPLEX){
		std::shared_ptr<MultiDouble> mp;
		try {
			switch(subSys){
			case 3:{ mp  = (paras.GetMultiDouble(0)); break; }
			case 4:{ mp  = (paras.GetMultiDouble(1)); break; }
			case 5:{ mp  = (paras.GetMultiDouble(2)); break; }
			}
		} catch ( std::exception &ex ) {
			BOOST_LOG_TRIVIAL(error) << "BreitWignerStrategy::execute() | "
					<<ex.what();
			throw( WrongParType("BreitWignerStrategy::execute() | "
					"Failed to obtain data vector from parameter list!")
			);
		}
		std::vector<std::complex<double> > results(mp->GetNValues(),
				std::complex<double>(0.));
		//calc BW for each point
		for(unsigned int ele=0; ele<mp->GetNValues(); ele++){
			double mSq = (mp->GetValue(ele));
			try{
				results[ele] = AmpRelBreitWignerRes::dynamicalFunction(
						mSq,m0,ma,mb,Gamma0,
						spin,d, formFactorType(ffType)
				);
			} catch ( std::exception& ex ) {
				BOOST_LOG_TRIVIAL(error) << "BreitWignerStrategy::execute() | "
						<<ex.what();
				throw( std::runtime_error("BreitWignerStrategy::execute() | "
						"Evaluation of dynamic function failed!")
				);
			}
		}
		out = std::shared_ptr<AbsParameter>(
				new MultiComplex(out->GetName(),results));
//		BOOST_LOG_TRIVIAL(debug) <<"BreitWignerStrategy::execute() | "
//				"finished!";
		return true;
	}//end multicomplex output

	//Only StandardDim Paras in input
	double mSq;
	try {
		switch(subSys){
		case 3:{ mSq  = paras.GetDoubleParameter(10)->GetValue(); break; }
		case 4:{ mSq  = paras.GetDoubleParameter(11)->GetValue(); break; }
		case 5:{ mSq  = paras.GetDoubleParameter(12)->GetValue(); break; }
		}
	} catch ( std::exception &ex ) {
		BOOST_LOG_TRIVIAL(error) << "BreitWignerStrategy::execute() | "
				<<ex.what();
		throw( WrongParType("BreitWignerStrategy::execute() | "
				"Failed to obtain data from parameter list!")
		);
	}
	std::complex<double> result;
	try{
		result = AmpRelBreitWignerRes::dynamicalFunction(
				mSq,m0,ma,mb,Gamma0, spin,d, formFactorType(ffType)
		);
	} catch (std::exception& ex) {
		BOOST_LOG_TRIVIAL(error) << "BreitWignerStrategy::execute() | "
				<<ex.what();
		throw( std::runtime_error("BreitWignerStrategy::execute() | "
				"Evaluation of dynamic function failed!")
		);
	}
	out = std::shared_ptr<AbsParameter>(
			new ComplexParameter(out->GetName(), result));
//	BOOST_LOG_TRIVIAL(debug) <<"BreitWignerStrategy::execute() | finished!";
	return true;
}
