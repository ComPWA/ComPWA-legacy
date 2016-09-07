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
		unsigned int varIdA, unsigned int varIdB,
		std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
		std::shared_ptr<DoubleParameter> mass,
		Spin spin, Spin m, Spin n, int P, int C,
		std::string mother, std::string particleA, std::string particleB,
		std::shared_ptr<DoubleParameter> width,
		std::shared_ptr<DoubleParameter> mesonRadius,
		std::shared_ptr<DoubleParameter> motherRadius,
		formFactorType type,
		int nCalls, normStyle nS) :
		AmpAbsDynamicalFunction(name, varIdA, varIdB, mag, phase, mass,
				spin, m, n, P, C,
				mother, particleA, particleB,
				mesonRadius, motherRadius, type, nCalls, nS),
				_width(width)
{
	if( _width->GetValue() != tmp_width) {
		SetModified();
		tmp_width = _width->GetValue();
	}
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

	return;
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
double AmpRelBreitWignerRes::GetIntegral()
{
	CheckModified();
	if(_modified){
		tmp_integral = integral();
		_modified=0;
	}
	return tmp_integral;
}

std::complex<double> AmpRelBreitWignerRes::EvaluateAmp(dataPoint& point)
{
	double mSq = point.getVal(_subSys);
	std::complex<double> result;
	try{
		result = dynamicalFunction(
				mSq,
				_mass->GetValue(),
				_mass1,
				_mass2,
				_width->GetValue(),
				_spin,
				_mesonRadius->GetValue(),
				_ffType
		);
	} catch (std::exception& ex){
		BOOST_LOG_TRIVIAL(error) <<"AmpRelBreitWignerRes::EvaluateAmp() | "
				"Dynamical function can not be evalutated: "<<ex.what();
		throw;
	}
	return result;
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
		ParameterList& sample, ParameterList& toySample,std::string suffix)
{
	DalitzKinematics* kin =
			dynamic_cast<DalitzKinematics*>(Kinematics::instance());
//	auto var1_limit = kin->GetMinMax( GetVarIdA() );
//	auto var2_limit = kin->GetMinMax( GetVarIdB() );
//	double phspVol = (var1_limit.second-var1_limit.first)
//			*(var2_limit.second-var2_limit.first);
	double phspVol = kin->GetPhspVolume();
//	double phspVol = 1;

	int sampleSize = sample.GetMultiDouble(0)->GetNValues();
	int toySampleSize = toySample.GetMultiDouble(0)->GetNValues();

	BOOST_LOG_TRIVIAL(info) << "AmpRelBreitWignerRes::setupBasicTree() | "
			<<_name << " nEvents=" <<sampleSize<<" nPhspEvents="
			<<toySampleSize;

	//------------Setup Tree---------------------
	std::shared_ptr<FunctionTree> newTree(new FunctionTree());

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
			new BreitWignerStrategy(_name) );
	std::shared_ptr<WignerDStrategy> angdStrat(
			new WignerDStrategy(_name) );

	//Reso=BW*C*AD*N
	newTree->createHead("Reso_"+_name, mmultStrat, sampleSize);

	newTree->createNode("PreFactor_"+_name, complStrat, "Reso_"+_name);
	newTree->createLeaf(
			"IntensPre_"+_name, std::abs(_prefactor), "PreFactor_"+_name);
	newTree->createLeaf(
			"PhasePre_"+_name, std::arg(_prefactor), "PreFactor_"+_name);

	newTree->createNode("C_"+_name, complStrat, "Reso_"+_name); //m0c
	newTree->createLeaf("Intens_"+_name, _mag, "C_"+_name); //r
	newTree->createLeaf("Phase_"+_name, _phase, "C_"+_name); //phi
	//Angular distribution
	if( _spin )
		newTree->insertTree(_wignerD.SetupTree(sample,suffix), "Reso_"+_name);

	//Breit-Wigner
	newTree->createNode("RelBW_"+_name, rbwStrat, "Reso_"+_name, sampleSize);
	newTree->createLeaf("mass", _mass, "RelBW_"+_name); //m0
	newTree->createLeaf("width", _width, "RelBW_"+_name); //resWidth
	newTree->createLeaf("spin", _spin, "RelBW_"+_name); //spin
	newTree->createLeaf("mesonRadius", _mesonRadius, "RelBW_"+_name); //d
	newTree->createLeaf("formFactorType", _ffType , "RelBW_"+_name); //d
	newTree->createLeaf("ma", _mass1, "RelBW_"+_name); //ma
	newTree->createLeaf("mb", _mass2, "RelBW_"+_name); //mb
	newTree->createLeaf("sample", sample.GetMultiDouble(_subSys), "RelBW_"+_name); //mc

	//Normalization
	if(_normStyle==normStyle::none) {
		newTree->createLeaf("N_"+_name, 1., "Reso_"+_name);
	} else {
		//Normalization parameter for dynamical amplitude
		newTree->createNode("N_"+_name, sqRootStrat, "Reso_"+_name); //N = sqrt(NSq)
		newTree->createNode("NSq_"+_name, multDStrat, "N_"+_name); //NSq = (PhspVolume/N_phspMC * Sum(|A|^2))^-1
		newTree->createLeaf("PhspSize_"+_name, toySampleSize, "NSq_"+_name); // N_phspMC
		newTree->createLeaf("PhspVolume_"+_name, 1/phspVol, "NSq_"+_name); // 1/PhspVolume
		newTree->createNode("InvSum_"+_name, invStrat, "NSq_"+_name); //1/Sum(|A|^2)
		newTree->createNode("Sum_"+_name, addStrat, "InvSum_"+_name); //Sum(|A|^2)
		newTree->createNode("AbsVal_"+_name, msqStrat, "Sum_"+_name); //|A_i|^2

		newTree->createNode("NormReso_"+_name, mmultStrat, "AbsVal_"+_name,
				toySampleSize);
		//Angular distribution (Normalization)
		if( _spin )
			newTree->insertTree(_wignerD.SetupTree(toySample,suffix),
					"NormReso_"+_name);
		//Breit-Wigner (Normalization)
		newTree->createNode("NormBW_"+_name, rbwStrat, "NormReso_"+_name,
				toySampleSize); //BW
		newTree->createLeaf("mass", _mass, "NormBW_"+_name); //m0
		newTree->createLeaf("width", _width, "NormBW_"+_name); //resWidth
		newTree->createLeaf("spin", _spin, "NormBW_"+_name); //spin
		newTree->createLeaf("mesonRadius", _mesonRadius, "NormBW_"+_name); //d
		newTree->createLeaf("formFactorType", _ffType , "NormBW_"+_name); //d
		newTree->createLeaf("ma", _mass1, "NormBW_"+_name); //ma
		newTree->createLeaf("mb", _mass2, "NormBW_"+_name); //mb
		newTree->createLeaf("phspSample", toySample.GetMultiDouble(_subSys),
				"NormBW_"+_name);
	}
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
	if( paras.GetNDouble() != 7 && paras.GetNDouble() != 8)
		throw( BadParameter("BreitWignerStrategy::execute() | "
				"number of DoubleParameters does not match!")
		);

	double Gamma0, m0, d, ma, mb;
	unsigned int spin;
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
	ma = paras.GetDoubleParameter(5)->GetValue();
	mb = paras.GetDoubleParameter(6)->GetValue();

	//	BOOST_LOG_TRIVIAL(debug) << "BreitWignerStrategy::execute() | mR="<<m0
	//			<<" Gamma="<<Gamma0<<" spin="<<spin<<" radius="<<d<<" ffType="<<ffType
	//			<<" subSys="<<subSys<<" ma="<<ma<<" mb="<<mb;

	//MultiDim output, must have multidim Paras in input
	if(checkType == ParType::MCOMPLEX){
		std::shared_ptr<MultiDouble> mp = paras.GetMultiDouble(0);
		std::vector<std::complex<double> > results(mp->GetNValues(),
				std::complex<double>(0.));
		//calc BW for each point
		for(unsigned int ele=0; ele<mp->GetNValues(); ele++){
			try{
				results.at(ele) = AmpRelBreitWignerRes::dynamicalFunction(
						mp->GetValue(ele),
						m0,ma,mb,Gamma0,
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
	double mSq = paras.GetDoubleParameter(7)->GetValue();
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
