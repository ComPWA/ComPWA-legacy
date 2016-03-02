//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------

#include <stdlib.h>
#include "gsl/gsl_monte.h"
#include "gsl/gsl_monte_plain.h"
#include "gsl/gsl_monte_miser.h"
#include "gsl/gsl_monte_vegas.h"

#include "Core/PhysConst.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Physics/AmplitudeSum/AmpAbsDynamicalFunction.hpp"

AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(const char *name,
		std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
		std::shared_ptr<DoubleParameter> mass, int part1, int part2,
		Spin spin, Spin m, Spin n,
		std::shared_ptr<DoubleParameter> mesonR, //  meson radius
		std::shared_ptr<DoubleParameter> motherR, //  mother radius
		formFactorType type,
		int nCalls, normStyle nS) :
		_name(name), _mag(mag), _phase(phase), _mass(mass), _subSys(part1+part2),
		_part1(part1), _part2(part2), _spin(spin),
		_m(m), _n(n), _mesonRadius(mesonR), _motherRadius(motherR), _ffType(type),
		_nCalls(nCalls), _normStyle(nS), _norm(1.0), modified(1),
		_wignerD(part1+part2, spin)
{
	initialize();
}

AmpAbsDynamicalFunction::AmpAbsDynamicalFunction(const char *name,
		std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
		std::shared_ptr<DoubleParameter> mass, int part1, int part2,
		Spin spin, Spin m, Spin n,
		formFactorType type,
		int nCalls, normStyle nS) :
		_name(name), _mag(mag), _phase(phase), _mass(mass), _part1(part1), _part2(part2),
		_subSys(part1+part2), _spin(spin), _m(m), _n(n),
		_mesonRadius(std::make_shared<DoubleParameter>(name, 1.0)),
		_motherRadius(std::make_shared<DoubleParameter>(name, 1.0)), _ffType(type),
		_nCalls(nCalls), _normStyle(nS), _norm(1.0), modified(1),
		_wignerD(part1+part2, spin)
{
	initialize();
}

std::string AmpAbsDynamicalFunction::to_str() const
{
	std::stringstream str;
	str<<"AmpAbsDynamicalFunction | "<<_name<<" enabled="<<_enable
			<< " nCalls="<<_nCalls<<" norm="<<_normStyle
			<< " subSys="<<_subSys<<" J="<<_spin<<" ffType="<<_ffType<<std::endl;
	str<<"Parameters:"<<std::endl;
	str<<_mag->to_str()<<std::endl;
	str<<_phase->to_str()<<std::endl;
	str<<_mass->to_str()<<std::endl;
	str<<_mesonRadius->to_str()<<std::endl;
	str<<_motherRadius->to_str()<<std::endl;

	return str.str();
}

void AmpAbsDynamicalFunction::Configure(
		boost::property_tree::ptree::value_type const& v,
		ParameterList& list)
{
	boost::property_tree::ptree pt = v.second;

	//Name (mandatory)
	auto tmp_name= pt.get_optional<std::string>("<xmlattr>.name");
	if(!tmp_name)
		tmp_name= pt.get_optional<std::string>("name");
	if(!tmp_name)
		throw BadParameter("AmpAbsDynamicalFunction::Configure() | "
				"No name for resonance specified!");
	_name = tmp_name.get();

	//Enable/Disable resonance (optional)
	auto tmp_enable= pt.get_optional<bool>("<xmlattr>.enable");
	if(!tmp_enable)
		tmp_enable= pt.get_optional<bool>("enable");
	if(!tmp_enable) _enable = 0;
	else _enable = tmp_enable.get();

	//Mother radius (optional)
	auto tmp_motherRadius_fix = pt.get<bool>("motherRadius_fix",1);
	auto tmp_motherRadius_min = pt.get<double>("motherRadius_min",0.);
	auto tmp_motherRadius_max = pt.get<double>("motherRadius_max",10.);
	auto tmp_motherRadius_name = pt.get_optional<std::string>("motherRadius_name");
	if(!tmp_motherRadius_name){
		//mother radius is not a strict requriement
		double tmp_motherRadius= pt.get<double>("motherRadius",1.0);
		_motherRadius = std::shared_ptr<DoubleParameter>(
				new DoubleParameter(
						"motherRadius_"+_name,tmp_motherRadius,
						tmp_motherRadius_min, tmp_motherRadius_max
				)
		);
		_motherRadius->FixParameter(tmp_motherRadius_fix);
		//if(_enable) list.AddParameter(_motherRadius);
	} else {
		try{
			_motherRadius = list.GetDoubleParameter(tmp_motherRadius_name.get());
		} catch (BadParameter& ex){
			BOOST_LOG_TRIVIAL(error) <<"AmpAbsDynamicalFunction::Configure() | "
					"Requesting parameter "<<tmp_motherRadius_name.get()<<" but"
							" was not found in parameter list. "
							"Continue since parameter is not mandatory!";
		}
			_motherRadius =	std::shared_ptr<DoubleParameter>(
							new DoubleParameter("motherRadius_"+_name,1.0)
			);
	}

	//Meson radius (mandatory)
	auto tmp_mesonRadius_fix = pt.get<bool>("mesonRadius_fix",1);
	auto tmp_mesonRadius_min = pt.get<double>("mesonRadius_min",0.0);
	auto tmp_mesonRadius_max = pt.get<double>("mesonRadius_max",10.0);
	auto tmp_mesonRadius_name = pt.get_optional<std::string>("mesonRadius_name");
	if(!tmp_mesonRadius_name){
		auto tmp_mesonRadius = pt.get_optional<double>("mesonRadius");
		if(!tmp_mesonRadius)
			throw BadParameter("AmpAbsDynamicalFunction::Configure() | "
					"mesonRadius for "+_name+" not specified!");
		_mesonRadius = std::shared_ptr<DoubleParameter>(
				new DoubleParameter(
						"d_"+_name,tmp_mesonRadius.get(),
						tmp_mesonRadius_min, tmp_mesonRadius_max
				)
		);
		_mesonRadius->FixParameter(tmp_mesonRadius_fix);
		if(_enable) list.AddParameter(_mesonRadius);
	} else {
		try{
			_mesonRadius = list.GetDoubleParameter(tmp_mesonRadius_name.get());
		} catch (BadParameter& ex){
			BOOST_LOG_TRIVIAL(error) <<"AmpAbsDynamicalFunction::Configure() | "
					"Requesting parameter "<<tmp_mesonRadius_name.get()<<" but"
							" was not found in parameter list. "
							"Quit since parameter is mandatory!";
			throw;
		}
	}

	//Mass (mandatory)
	auto tmp_mass_fix = pt.get<bool>("mass_fix",1);
	auto tmp_mass_min = pt.get<double>("mass_min",0.0);
	auto tmp_mass_max = pt.get<double>("mass_max",10.0);
	auto tmp_mass_name = pt.get_optional<std::string>("mass_name");
	if(!tmp_mass_name){
		auto tmp_mass = pt.get_optional<double>("mass");
		if(!tmp_mass)
			throw BadParameter("AmpAbsDynamicalFunction::Configure() | "
					"mass for "+_name+" not specified!");
		_mass = std::shared_ptr<DoubleParameter>(
				new DoubleParameter(
						"mass_"+_name,tmp_mass.get(),
						tmp_mass_min, tmp_mass_max
				)
		);
		_mass->FixParameter(tmp_mass_fix);
		if(_enable) list.AddParameter(_mass);
	} else {
		try{
			_mass = list.GetDoubleParameter(tmp_mass_name.get());
		} catch (BadParameter& ex){
			BOOST_LOG_TRIVIAL(error) <<"AmpAbsDynamicalFunction::Configure() | "
					"Requesting parameter "<<tmp_mass_name.get()<<" but"
							" was not found in parameter list. "
							"Quit since parameter is mandatory!";
			throw;
		}
	}

	//Strength (mandatory)
	auto tmp_strength_fix = pt.get<bool>("strength_fix",1);
	auto tmp_strength_min = pt.get<double>("strength_min",0.0);
	auto tmp_strength_max = pt.get<double>("strength_max",5.0);
	auto tmp_strength_name = pt.get_optional<std::string>("strength_name");
	if(!tmp_strength_name){
		auto tmp_strength = pt.get_optional<double>("strength");
		if(!tmp_strength)
			throw BadParameter("AmpAbsDynamicalFunction::Configure() | "
					"strength for "+_name+" not specified!");
		_mag = std::shared_ptr<DoubleParameter>(
				new DoubleParameter(
						"mag_"+_name,tmp_strength.get(),
						tmp_strength_min, tmp_strength_max
				)
		);
		_mag->FixParameter(tmp_strength_fix);
		if(_enable) list.AddParameter(_mag);
	} else {
		try{
			_mag = list.GetDoubleParameter(tmp_strength_name.get());
		} catch (BadParameter& ex){
			BOOST_LOG_TRIVIAL(error) <<"AmpAbsDynamicalFunction::Configure() | "
					"Requesting parameter "<<tmp_strength_name.get()<<" but"
							" was not found in parameter list. "
							"Quit since parameter is mandatory!";
			throw;
		}
	}

	//Phase (mandatory)
	auto tmp_phase_fix = pt.get<bool>("phase_fix",1);
	auto tmp_phase_min = pt.get<double>("phase_min",-300.0);
	auto tmp_phase_max = pt.get<double>("phase_max",300.0);
	auto tmp_phase_name = pt.get_optional<std::string>("phase_name");
	if(!tmp_phase_name){
		auto tmp_phase= pt.get_optional<double>("phase");
		if(!tmp_phase)
			throw BadParameter("AmpAbsDynamicalFunction::Configure() | "
					"phase for "+_name+" not specified!");
		_phase = std::shared_ptr<DoubleParameter>(
				new DoubleParameter(
						"phase_"+_name,tmp_phase.get(),
						tmp_phase_min, tmp_phase_max
				)
		);
		_phase->FixParameter(tmp_phase_fix);
		if(_enable) list.AddParameter(_phase);
	} else {
		try{
			_phase = list.GetDoubleParameter(tmp_phase_name.get());
		} catch (BadParameter& ex){
			BOOST_LOG_TRIVIAL(error) <<"AmpAbsDynamicalFunction::Configure() | "
					"Requesting parameter "<<tmp_phase_name.get()<<" but"
							" was not found in parameter list. "
							"Quit since parameter is mandatory!";
			throw;
		}
	}

	//FormFactor (optional)
	auto tmp_ffType= pt.get<unsigned int>("FormFactorType", 1);
	_ffType = formFactorType(tmp_ffType);

	auto tmp_spin = pt.get_optional<int>("spin");
	if(!tmp_spin)
		throw BadParameter("AmpAbsDynamicalFunction::Configure() | "
				"spin for "+_name+" not specified!");
	_spin = Spin(tmp_spin.get());

	auto tmp_part1 = pt.get_optional<int>("daughterA");
	if(!tmp_part1)
		throw BadParameter("AmpAbsDynamicalFunction::Configure() | "
				"daughterA for "+_name+" not specified!");
	_part1 = tmp_part1.get();
	auto tmp_part2 = pt.get_optional<int>("daughterB");
	if(!tmp_part2)
		throw BadParameter("AmpAbsDynamicalFunction::Configure() | "
				"daughterB for "+_name+" not specified!");
	_part2 = tmp_part2.get();

	_subSys = _part1+_part2;

	//optional parameters
	double tmp_m = pt.get<int>("m",0);
	_m = Spin(tmp_m);
	double tmp_n = pt.get<int>("n",0);
	_n = Spin(tmp_n);



	_wignerD = AmpWigner2(_subSys,_spin);
	initialize();
}

void AmpAbsDynamicalFunction::put(boost::property_tree::ptree &pt){
	pt.put("<xmlattr>.name", _name);
	pt.put("<xmlattr>.enable", _enable);
	pt.put("strength", _mag->GetValue());
	pt.put("strength_fix", _mag->IsFixed());
	pt.put("strength_min", _mag->GetMinValue());
	pt.put("strength_max", _mag->GetMaxValue());
	pt.put("phase", _phase->GetValue());
	pt.put("phase_fix", _phase->IsFixed());
	pt.put("phase_min", _phase->GetMinValue());
	pt.put("phase_max", _phase->GetMaxValue());
	pt.put("mass", _mass->GetValue());
	pt.put("mass_fix", _mass->IsFixed());
	pt.put("mass_min", _mass->GetMinValue());
	pt.put("mass_max", _mass->GetMaxValue());
	pt.put("FormFactorType", _ffType);
	pt.put("daughterA", _part1);
	pt.put("daughterB", _part2);
}

void AmpAbsDynamicalFunction::CheckModified() {
	if(_mass->GetValue() != tmp_mass){
		SetModified();
		tmp_mass = _mass->GetValue();
	}
}

double AmpAbsDynamicalFunction::GetInvMass(dataPoint& point){
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double mSq = -999;
	switch(_subSys){
	case 3: mSq=kin->getThirdVariableSq(point.getVal(0),point.getVal(1)); break;
	case 4: mSq=point.getVal(1); break;
	case 5: mSq=point.getVal(0); break;
	}
	return mSq;
}

void AmpAbsDynamicalFunction::initialize()
{
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	_M=kin->M;
	if(_subSys==5){
		_ma=kin->m3;
		_mb=kin->m2;
	}
	if(_subSys==4){
		_ma=kin->m3;
		_mb=kin->m1;
	}
	if(_subSys==3){
		_ma=kin->m2;
		_mb=kin->m1;
	}

	if(_mass->GetValue() != tmp_mass){
		SetModified();
		tmp_mass = _mass->GetValue();
	}
}

AmpAbsDynamicalFunction::~AmpAbsDynamicalFunction() 
{
}

std::complex<double> AmpAbsDynamicalFunction::GetCoefficient() {
	return std::complex<double>(
			std::fabs(_mag->GetValue())*cos(_phase->GetValue()),
			std::fabs(_mag->GetValue())*sin(_phase->GetValue())
	);
}

std::complex<double> AmpAbsDynamicalFunction::evaluate(dataPoint& point){
	CheckModified();
	std::complex<double> res = evaluateAmp(point);
	double ang = evaluateWignerD(point);
	return (GetCoefficient()*GetNormalization()*res*ang);
}

double evalAmp(double* x, size_t dim, void* param) {
	/* We need a wrapper here because a eval() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	dataPoint pp; pp.setVal(0,x[1]);pp.setVal(1,x[0]);
	if( !kin->isWithinPhsp(pp) ) return 0;//only integrate over phase space
	std::complex<double> res = static_cast<AmpAbsDynamicalFunction*>(param)->evaluateAmp(pp);
	//include angular distribution in normalization
	res *= static_cast<AmpAbsDynamicalFunction*>(param)->evaluateWignerD(pp);
	return ( std::norm(res) ); //integrate over |F|^2
}

double AmpAbsDynamicalFunction::integral(){
	size_t dim=2;
	double res=0.0, err=0.0;

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function F = {&evalAmp,dim, const_cast<AmpAbsDynamicalFunction*> (this)};
	//	gsl_monte_function F = {&twoDimGaussian,dim, new int()};//using test function; result should be 1

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&F, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	BOOST_LOG_TRIVIAL(debug)<<"AmpAbsDynamicalFunction::integral() Integration result for |"
			<<_name<<"|^2: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

double AmpAbsDynamicalFunction::GetNormalization(){
	if(_norm<0) return 1.0; //normalization is disabled
	//	return _norm; //disable recalculation of normalization
	if(!modified) return _norm;
	_norm = 1/sqrt(integral());
	modified=0;
	return _norm;
}

double eval(double* x, size_t dim, void* param) {
	/* We need a wrapper here because evaluate() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	dataPoint pp; pp.setVal(0,x[1]);pp.setVal(1,x[0]);
	if( !kin->isWithinPhsp(pp) ) return 0;//only integrate over phase space
	std::complex<double> res = static_cast<AmpAbsDynamicalFunction*>(param)->evaluateAmp(pp);
	double ang = static_cast<AmpAbsDynamicalFunction*>(param)->evaluateWignerD(pp);
	double norm = static_cast<AmpAbsDynamicalFunction*>(param)->GetNormalization();
	return ( std::norm(res*ang*norm) ); //integrate over |F|^2
}

double AmpAbsDynamicalFunction::totalIntegral() const{
	//Save CPU time
	return 1;

	size_t dim=2;
	double res=0.0, err=0.0;

	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function F = {&eval,dim, const_cast<AmpAbsDynamicalFunction*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&F, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	//	BOOST_LOG_TRIVIAL(debug)<<"AmpAbsDynamicalFunction::totalIntegral() result for |"<<_name<<"|^2: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

std::complex<double> AmpAbsDynamicalFunction::widthToCoupling(double mSq, double mR, double width,
		double ma, double mb, double spin, double mesonRadius, formFactorType type)
{
	double sqrtS = sqrt(mSq);
	//calculate gammaA(s_R)
	double ffR = Kinematics::FormFactor(mR,ma,mb,spin,mesonRadius,type);
	std::complex<double> qR = Kinematics::qValue(mR,ma,mb);
	//calculate phsp factor
	std::complex<double> rho = Kinematics::phspFactor(sqrtS,ma,mb);
	std::complex<double> denom = std::pow(qR,spin)*ffR*sqrt(rho);
	std::complex<double> result = std::complex<double>(sqrt(mR*width), 0) / denom;
	return result;
}

std::complex<double> AmpAbsDynamicalFunction::couplingToWidth(double mSq, double mR, double g,
		double ma, double mb, double spin, double mesonRadius, formFactorType type)
{
	double sqrtM = sqrt(mSq);
	//calculate gammaA(s_R)
	double ffR = Kinematics::FormFactor(mR,ma,mb,spin,mesonRadius,type);
	std::complex<double> qR = std::pow(Kinematics::qValue(mR,ma,mb),spin);
	std::complex<double> gammaA = ffR*qR;
	//calculate phsp factor
	std::complex<double> rho = Kinematics::phspFactor(sqrtM,ma,mb);
	std::complex<double> result = std::norm(gammaA)*g*g*rho/ mR;
	if(result.real()!=result.real() || result.imag()!=result.imag()){
		std::cout<<"AmpKinematics::couplingToWidth() | NaN! mSq="<<mSq
				<<" mR="<<mR<<" g="<<g<<" ma="<<ma<<" mb="<<mb<<std::endl;
		std::cout<<qR<<" "<<gammaA<<" "<<rho<<" "<<g<<std::endl;
	}
	return result;
}
double twoDimGaussian(double* z, size_t dim, void *param){
	if(dim!=2) return 0;
	/* test environment for numeric integration:
	 * 	Calculating integral of normalized gaussian:
	 * 	f(x,y) = A exp( - (x-x0)^2/(2*sigmaX^2) + (y-y0)^2/(2*sigmaY^2)
	 * 	with A=1/(2*pi*sigmaX*sigmaY) this function is normalized to 1
	 */
	double x = z[0]; double y = z[1];
	//mean and width need to be adjusted according to final state kinematics
	double x0=1.1, y0=1.1; //mean
	double sigmaX=0.01, sigmaY=0.01; //width
	double pi = PhysConst::instance()->getConstValue("Pi");

	double result = exp( -(x-x0)*(x-x0)/(2*sigmaX*sigmaX) - (y-y0)*(y-y0)/(2*sigmaY*sigmaY) );
	result/=2*pi*sigmaY*sigmaX;
	return result;
}
