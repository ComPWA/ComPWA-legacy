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

#ifndef AMP_REL_BREIT_WIGNER_RES
#define AMP_REL_BREIT_WIGNER_RES

#include <vector>

#include "Core/Functions.hpp"
#include "Core/Exceptions.hpp"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"
#include "Physics/AmplitudeSum/AmpKinematics.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"
#include "Physics/AmplitudeSum/NonResonant.hpp"

class AmpRelBreitWignerRes : public AmpAbsDynamicalFunction, public AmpKinematics {
public:

	AmpRelBreitWignerRes(const char *name,
			std::shared_ptr<DoubleParameter> _resMass, std::shared_ptr<DoubleParameter> _resWidth,
			std::shared_ptr<DoubleParameter> _radius, std::shared_ptr<DoubleParameter> _motherRadius,
			int _subsys, int resSpin, int m, int n
	) ;
	AmpRelBreitWignerRes(const AmpRelBreitWignerRes&, const char*);
	AmpRelBreitWignerRes(const AmpRelBreitWignerRes&);

	virtual ~AmpRelBreitWignerRes();

	static std::complex<double> dynamicalFunction(double mSq, double mR, double ma, double mb, double gamma0, unsigned int J, double mesonRadius);
	virtual void initialise();
	virtual std::complex<double> evaluate(dataPoint& point) { return _norm*evaluateAmp(point)*evaluateWignerD(point); };
	virtual std::complex<double> evaluateAmp(dataPoint& point);
//	virtual double evaluateWignerD(dataPoint& point) {
//		if(_spin==0) return 1.0;//save some computing time
//		return _wignerD.evaluate(point);
//	};

	void setDecayMasses(double, double, double, double);
	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	unsigned int getNParams(){ return nParams;}

protected:
	std::shared_ptr<DoubleParameter> _resWidth;
//	AmpWigner2 _wignerD;
	bool foundMasses;
	unsigned int id23, id13;
	unsigned int nParams;

};

class BreitWignerConf : public basicConf
{
public:
	BreitWignerConf(const boost::property_tree::ptree &pt_) : basicConf(pt_){
		m_mass= pt_.get<double>("mass");
		m_mass_fix= pt_.get<bool>("mass_fix");
		m_mass_min= pt_.get<double>("mass_min");
		m_mass_max= pt_.get<double>("mass_max");
		m_width= pt_.get<double>("width");
		m_width_fix= pt_.get<bool>("width_fix");
		m_width_min= pt_.get<double>("width_min");
		m_width_max= pt_.get<double>("width_max");
		m_mesonRadius= pt_.get<double>("mesonRadius");
		m_spin= pt_.get<unsigned int>("spin");
		m_m= pt_.get<unsigned int>("m");
		m_n= pt_.get<unsigned int>("n");
		m_daughterA= pt_.get<unsigned int>("daughterA");
		m_daughterB= pt_.get<unsigned int>("daughterB");
	}
	virtual void put(boost::property_tree::ptree &pt_){
		basicConf::put(pt_);
		pt_.put("mass", m_mass);
		pt_.put("mass_fix", m_mass_fix);
		pt_.put("mass_min", m_mass_min);
		pt_.put("mass_max", m_mass_max);
		pt_.put("width", m_width);
		pt_.put("width_fix", m_width_fix);
		pt_.put("width_min", m_width_min);
		pt_.put("width_max", m_width_max);
		pt_.put("mesonRadius", m_mesonRadius);
		pt_.put("spin", m_spin);
		pt_.put("m", m_m);
		pt_.put("n", m_n);
		pt_.put("daughterA", m_daughterA);
		pt_.put("daughterB", m_daughterB);
	}
	virtual void update(ParameterList par){
		basicConf::update(par);
		m_mass= par.GetDoubleParameter("m0_"+m_name)->GetValue();
		m_width= par.GetDoubleParameter("width_"+m_name)->GetValue();
	}

	double m_mass;
	bool m_mass_fix;
	double m_mass_min;
	double m_mass_max;
	double m_width;
	bool m_width_fix;
	double m_width_min;
	double m_width_max;

	double m_mesonRadius;
	unsigned int m_spin;
	unsigned int m_m;
	unsigned int m_n;

	unsigned int m_daughterA;
	unsigned int m_daughterB;
};

class BreitWignerStrategy : public Strategy {
public:
	BreitWignerStrategy(const std::string resonanceName, ParType in):Strategy(in),name(resonanceName){}

	virtual const std::string to_str() const {
		return ("relativistic BreitWigner of "+name);
	}

	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
		if( checkType != out->type() ) {
			throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
			return false;
		}

		double Gamma0, m0, d, ma, mb;
		unsigned int spin, subSys;
		//Get parameters from ParameterList -
		//enclosing in try...catch for the case that names of nodes have changed
		try{
			m0 = double(paras.GetParameterValue("m0_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerStrategy: can't find parameter m0_"+name;
			throw;
		}
		try{
			spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerStrategy: can't find parameter ParOfNode_spin_"+name;
			throw;
		}
		try{
			d = double(paras.GetParameterValue("d_"+name));
		}catch(BadParameter& e){
			std::cout<<paras<<std::endl;
			BOOST_LOG_TRIVIAL(error) <<"----BreitWignerStrategy: can't find parameter d_"+name;
			throw;
		}
		try{
			subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerStrategy: can't find parameter ParOfNode_subSysFlag_"+name;
			throw;
		}
		try{
			ma = double(paras.GetParameterValue("ParOfNode_ma_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerStrategy: can't find parameter ParOfNode_ma_"+name;
			throw;
		}
		try{
			mb = double(paras.GetParameterValue("ParOfNode_mb_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerStrategy: can't find parameter ParOfNode_mb_"+name;
			throw;
		}

		try{
			Gamma0 = double(paras.GetParameterValue("width_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerStrategy: can't find parameter width_"+name;
			throw;
		}

		//MultiDim output, must have multidim Paras in input
		if(checkType == ParType::MCOMPLEX){
			if(paras.GetNMultiDouble()){
				unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();

				std::vector<std::complex<double> > results(nElements, std::complex<double>(0.));
				std::shared_ptr<MultiDouble> mp;//=paras.GetMultiDouble("mym_"+name);
				switch(subSys){
				case 3:{ //reso in sys of particles 1&2
					mp  = (paras.GetMultiDouble("m12sq"));
					break;
				}
				case 4:{ //reso in sys of particles 1&3
					mp  = (paras.GetMultiDouble("m13sq"));
					break;
				}
				case 5:{ //reso in sys of particles 2&3
					mp  = (paras.GetMultiDouble("m23sq"));
					break;
				}
				}

				//calc BW for each point
				for(unsigned int ele=0; ele<nElements; ele++){
					double mSq = (mp->GetValue(ele));
					results[ele] = AmpRelBreitWignerRes::dynamicalFunction(mSq,m0,ma,mb,Gamma0,spin,d);
				}
				out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));
				return true;
			}else{ //end multidim para treatment
				throw(WrongParType("Input MultiDoubles missing in BW strat of "+name));
				return false;
			}
		}//end multicomplex output

		//Only StandardDim Paras in input
		//  double spinTerm = evaluateWignerD(); //spinTerm =1;
		double mSq;// = sqrt(paras.GetParameterValue("mym_"+name));
		switch(subSys){
		case 3:{ //reso in sys of particles 1&2
			mSq  = (double(paras.GetParameterValue("m12sq")));
			break;
		}
		case 4:{ //reso in sys of particles 1&3
			mSq  = (double(paras.GetParameterValue("m13sq")));
			break;
		}
		case 5:{ //reso in sys of particles 2&3
			mSq  = (double(paras.GetParameterValue("m23sq")));
			break;
		}
		}
		std::complex<double> result = AmpRelBreitWignerRes::dynamicalFunction(mSq,m0,ma,mb,Gamma0,spin,d);
		out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));
		return true;
	}

protected:
	std::string name;
};

class BreitWignerPhspStrategy : public BreitWignerStrategy {
public:
	BreitWignerPhspStrategy(const std::string resonanceName, ParType in):BreitWignerStrategy(resonanceName,in){}
	virtual bool execute(ParameterList& paras, std::shared_ptr<AbsParameter>& out) {
		if( checkType != out->type() ) {
			throw(WrongParType(std::string("Output Type ")+ParNames[out->type()]+std::string(" conflicts expected type ")+ParNames[checkType]+std::string(" of ")+name+" BW strat"));
			return false;
		}

		double Gamma0, m0, d, ma, mb;
		unsigned int spin, subSys;
		//Get parameters from ParameterList -
		//enclosing in try...catch for the case that names of nodes have changed
		try{
			m0 = double(paras.GetParameterValue("m0_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerPhspStrategy: can't find parameter m0_"+name;
			throw;
		}
		try{
		spin = (unsigned int)(paras.GetParameterValue("ParOfNode_spin_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerPhspStrategy: can't find parameter ParOfNode_spin_"+name;
			throw;
		}
		try{
		d = double(paras.GetParameterValue("d_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerPhspStrategy: can't find parameter d_"+name;
			throw;
		}
		try{
		subSys = double(paras.GetParameterValue("ParOfNode_subSysFlag_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerPhspStrategy: can't find parameter ParOfNode_subSysFlag_"+name;
			throw;
		}
		try{
		ma = double(paras.GetParameterValue("ParOfNode_ma_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerPhspStrategy: can't find parameter ParOfNode_ma_"+name;
			throw;
		}
		try{
		mb = double(paras.GetParameterValue("ParOfNode_mb_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerPhspStrategy: can't find parameter ParOfNode_mb_"+name;
			throw;
		}
		try{
			Gamma0 = double(paras.GetParameterValue("width_"+name));
		}catch(BadParameter& e){
			BOOST_LOG_TRIVIAL(error) <<"BreitWignerPhspStrategy: can't find parameter width_"+name;
			throw;
		}

		//MultiDim output, must have multidim Paras in input
		if(checkType == ParType::MCOMPLEX){
			if(paras.GetNMultiDouble()){
				unsigned int nElements = paras.GetMultiDouble(0)->GetNValues();

				std::vector<std::complex<double> > results(nElements, std::complex<double>(0.));
				std::shared_ptr<MultiDouble> mp;//=paras.GetMultiDouble("mym_"+name);
				switch(subSys){
				case 3:{ //reso in sys of particles 1&2
					mp  = (paras.GetMultiDouble("m12sq_phsp"));
					//map  = (paras.GetMultiDouble("m23"));
					// mbp  = (paras.GetMultiDouble("m13"));
					break;
				}
				case 4:{ //reso in sys of particles 1&3
					mp  = (paras.GetMultiDouble("m13sq_phsp"));
					// map  = (paras.GetMultiDouble("m12"));
					// mbp  = (paras.GetMultiDouble("m23"));
					break;
				}
				case 5:{ //reso in sys of particles 2&3
					mp  = (paras.GetMultiDouble("m23sq_phsp"));
					// map  = (paras.GetMultiDouble("m13"));
					// mbp  = (paras.GetMultiDouble("m12"));
					break;
				}
				}

				//calc BW for each point
				for(unsigned int ele=0; ele<nElements; ele++){
					double mSq = (mp->GetValue(ele));
					results[ele] = AmpRelBreitWignerRes::dynamicalFunction(mSq,m0,ma,mb,Gamma0,spin,d);
				}
				//std::vector<std::complex<double> > resultsTMP(nElements, std::complex<double>(1.));
				out = std::shared_ptr<AbsParameter>(new MultiComplex(out->GetName(),results));
				return true;
			}else{ //end multidim para treatment
				throw(WrongParType("Input MultiDoubles missing in BW strat of "+name));
				return false;
			}
		}//end multicomplex output


		//Only StandardDim Paras in input
		//  double spinTerm = evaluateWignerD(); //spinTerm =1;
		double mSq;// = sqrt(paras.GetParameterValue("mym_"+name));
		switch(subSys){
		case 3:{ //reso in sys of particles 1&2
			mSq  = (double(paras.GetParameterValue("m12sq_phsp")));
			break;
		}
		case 4:{ //reso in sys of particles 1&3
			mSq  = (double(paras.GetParameterValue("m13sq_phsp")));
			break;
		}
		case 5:{ //reso in sys of particles 2&3
			mSq  = (double(paras.GetParameterValue("m23sq_phsp")));
			break;
		}
		}
		std::complex<double> result = AmpRelBreitWignerRes::dynamicalFunction(mSq,m0,ma,mb,Gamma0,spin,d);
		std::cout<<"Strategy BWrelPHSP single value "<<result<<std::endl;

		out = std::shared_ptr<AbsParameter>(new ComplexParameter(out->GetName(), result));

		return true;
	}
}; // end BreitwignerPhspStrategy


#endif
