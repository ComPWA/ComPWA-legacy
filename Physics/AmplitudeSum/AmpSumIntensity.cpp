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

#include "Core/PhysConst.hpp"
#include "Core/Functions.hpp"

// Boost header files go here
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Physics/AmplitudeSum/NonResonant.hpp"
#include "Physics/AmplitudeSum/AmpRelBreitWignerRes.hpp"
#include "Physics/AmplitudeSum/AmpGausRes.hpp"
#include "Physics/AmplitudeSum/AmpFlatteRes.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

#include "gsl/gsl_monte_vegas.h"

#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"

AmpSumIntensity::AmpSumIntensity(normStyle ns, std::shared_ptr<Efficiency> eff,
		unsigned int nCalls) :
		_normStyle(ns), _calcMaxFcnVal(0), eff_(eff), _nCalls(nCalls)
{
	result.AddParameter(
			std::shared_ptr<DoubleParameter>(
					new DoubleParameter("AmpSumResult")
			)
	);
	return;
}

void AmpSumIntensity::LoadParameters(std::string path){
//	std::ifstream ifs3("ParList.xml");
//	boost::archive::xml_iarchive ia3(ifs3);
//	ia3 >> boost::serialization::make_nvp("ParList", list2);

}

//! Configure resonance from ptree
void AmpSumIntensity::Configure(const boost::property_tree::ptree &pt)
{
	BOOST_FOREACH( ptree::value_type const& v,
			pt.get_child("amplitude_setup") ) {
		/* We try to configure each type of resonance. In case that v.first does
		 * not contain the correct string, a BadConfig is thrown. We catch it
		 * and try the next type.
		 */
		try{
			auto amp=new AmpRelBreitWignerRes(_normStyle, _nCalls);
			amp->Configure(v,params);
			resoList.push_back(std::shared_ptr<AmpAbsDynamicalFunction>(amp));
			BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::Configure() | "
					"adding amplitude:"
					<<std::endl<<amp->to_str();
		} catch (BadConfig& ex) {}
		try{
			auto amp=new AmpFlatteRes(_normStyle, _nCalls);
			amp->Configure(v,params);
			resoList.push_back(std::shared_ptr<AmpAbsDynamicalFunction>(amp));
			BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::Configure() | "
					"adding amplitude:"
					<<std::endl<<amp->to_str();
		} catch (BadConfig& ex) {}
		try{
			auto amp=new NonResonant(_normStyle, _nCalls);
			amp->Configure(v,params);
			resoList.push_back(std::shared_ptr<AmpAbsDynamicalFunction>(amp));
			BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::Configure() | "
					"adding amplitude:"
					<<std::endl<<amp->to_str();
		} catch (BadConfig& ex) {}
	}
	BOOST_LOG_TRIVIAL(info) << "AmpSumIntensity::Configure() | "
			"Loaded property tree!";
	return;
}

//! Save resonance from to ptree
void AmpSumIntensity::Save(std::string filePath)
{
	using boost::property_tree::ptree;
	boost::property_tree::ptree ptSub, pt;

	BOOST_LOG_TRIVIAL(info) << "AmpSumIntensity::Save() | "
			"Save amplitude to "<<filePath<<"!";
	auto it = resoList.begin();
	for( ; it!= resoList.end(); ++it)
		(*it)->Save(ptSub);

	pt.add_child("amplitude_setup",ptSub);
	//new line at the end
	boost::property_tree::xml_writer_settings<char> settings('\t', 1);
	// Write the property tree to the XML file.
	write_xml(filePath, pt,std::locale(), settings);

}

std::shared_ptr<FunctionTree> AmpSumIntensity::setupBasicTree(
		allMasses& theMasses,allMasses& toyPhspSample,std::string suffix)
{
	BOOST_LOG_TRIVIAL(debug) << "AmpSumIntensity::setupBasicTree() generating new tree!";
	if(theMasses.nEvents==0){
		BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::setupBasicTree() data sample empty!";
		return std::shared_ptr<FunctionTree>();
	}
	if(toyPhspSample.nEvents==0){
		BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::setupBasicTree() toy sample empty!";
		return std::shared_ptr<FunctionTree>();
	}

	//------------Setup Tree---------------------
	std::shared_ptr<FunctionTree> newTree(new FunctionTree());

	//----Strategies needed
	std::shared_ptr<AddAll> maddStrat(new AddAll(ParType::MCOMPLEX));
	newTree->createHead("Amplitude"+suffix, maddStrat, theMasses.nEvents);

	auto it = resoList.begin();
	for( ; it!=resoList.end(); ++it){
		if(!(*it)->GetEnable()) continue;
		std::shared_ptr<FunctionTree> resTree= (*it)->setupTree(theMasses,
				toyPhspSample, ""+(*it)->GetName(),params);
		if(!resTree->sanityCheck())
			throw std::runtime_error("AmpSumIntensity::setupBasicTree() | "
					"Resonance tree didn't pass sanity check!");
		resTree->recalculate();
		newTree->insertTree(resTree, "Amplitude"+suffix);
	}

	BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::setupBasicTree(): tree constructed!!";
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

double AmpSumIntensity::evaluate(double x[], size_t dim) {
	/* Calculation amplitude integral (excluding efficiency) */
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
}

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
	BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::integrate() Integration result for amplitude sum: "
			<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

double interferenceIntegralWrapper(double* x, size_t dim, void* param) {
	/* Calculation amplitude integral (excluding efficiency) */
	if(dim!=2) return 0;
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	dataPoint point; point.setVal(1,x[0]); point.setVal(0,x[1]);
	if( !kin->isWithinPhsp(point) ) return 0;//only integrate over phase space

	AmpSumIntensity* amp = static_cast<AmpSumIntensity*>(param);
	std::shared_ptr<AmpAbsDynamicalFunction> amp1 = amp->GetResonance(amp->getFirstRes());
	std::shared_ptr<AmpAbsDynamicalFunction> amp2 = amp->GetResonance(amp->getSecondRes());
	std::complex<double> val1 = amp1->evaluate(point);
	std::complex<double> val2 = amp2->evaluate(point);

	//	if(amp->getFirstRes() != amp->getSecondRes())
	//		std::cout<<(val1*std::conj(val2)).real()
	//		<<" "<<(val2*std::conj(val1)).real()
	//		<< " "<<std::norm(val1)
	//		<< " "<<std::norm(val2)<<std::endl;
	double intens = (val1*std::conj(val2)).real();
	if(amp->getFirstRes() == amp->getSecondRes()) return intens;
	return 2*intens;
}

const double AmpSumIntensity::interferenceIntegral(unsigned int a, unsigned int b){
	assert( GetNumberOfResonances() > a );
	assert( GetNumberOfResonances() > b );
	firstRes = a; secondRes = b;

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
	gsl_monte_function G = {&interferenceIntegralWrapper,dim, const_cast<AmpSumIntensity*> (this)};

	/*	Choosing vegas algorithm here, because it is the most accurate:
	 * 		-> 10^5 calls gives (in my example) an accuracy of 0.03%
	 * 		 this should be sufficiency for most applications
	 */
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&G, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	BOOST_LOG_TRIVIAL(debug)<<"AmpSumIntensity::interferenceIntegrate() Interference term "
			"of "<<a<<" and "<<b<<" : "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

double AmpSumIntensity::evaluateEff(double x[], size_t dim) {
	/* Calculation amplitude normalization (including efficiency) */
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
	return normalization();
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
	//	BOOST_LOG_TRIVIAL(info)<<"AmpSumIntensity::normalization() Integration result for amplitude sum: "
	//			<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}

const double AmpSumIntensity::sliceIntensity(dataPoint& dataP, ParameterList& par,std::complex<double>* reso, unsigned int nResos){
	setParameterList(par);

	double AMPpdf=0;
	//TODO: implement slice fit
	//	if(Kinematics::instance()->isWithinPhsp(dataP)) AMPpdf = totAmp.evaluateSlice(dataP, reso, nResos,5);
	if(AMPpdf!=AMPpdf){
		BOOST_LOG_TRIVIAL(error)<<"Error AmpSumIntensity: Intensity is not a number!!";
		AMPpdf = 0;
	}
	double eff=eff_->evaluate(dataP);
	return AMPpdf*eff;
}

const ParameterList& AmpSumIntensity::intensity(std::vector<double> point, ParameterList& par){
	setParameterList(par);
	dataPoint dataP(point);
	return intensity(dataP);
}

const ParameterList& AmpSumIntensity::intensity(dataPoint& point, ParameterList& par){
	setParameterList(par);
	return intensity(point);
}

const ParameterList& AmpSumIntensity::intensityNoEff(dataPoint& point){
	std::complex<double> AMPpdf(0,0);
	if(Kinematics::instance()->isWithinPhsp(point)) {
		auto it = resoList.begin();
		for( ; it!=resoList.end(); ++it){
			if(!(*it)->GetEnable()) continue;
			AMPpdf += (*it)->evaluate(point);
		}
	}
	result.SetParameterValue(0,std::norm(AMPpdf));
	return result;
}

const ParameterList& AmpSumIntensity::intensity(dataPoint& point){
	intensityNoEff(point);
	double ampNoEff = result.GetParameterValue(0);
	double eff=eff_->evaluate(point);
	result.SetParameterValue(0,ampNoEff*eff);
	return result;
}

void AmpSumIntensity::setParameterList(ParameterList& par){
	//parameters varied by Minimization algorithm
	if(par.GetNDouble()!=params.GetNDouble())
		throw std::runtime_error("AmpSumIntensity::setParameterList(): size of parameter lists don't match");
	//Should we compared the parameter names? String comparison is slow
	for(unsigned int i=0; i<params.GetNDouble(); i++)
		params.GetDoubleParameter(i)->UpdateParameter(par.GetDoubleParameter(i));
	return;
}

bool AmpSumIntensity::copyParameterList(ParameterList& outPar){
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
			outStr<<"-------- "<<GetNameOfResonance(n)<<" ---------\n";
			n++;
		}
		outStr<<p->GetName()<<" = "<<p->GetValue();
		if(p->HasError())
			outStr<<"+-"<<p->GetError();
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
	for(unsigned int i=0;i<resoList.size();i++){
		if(!resoList.at(i)->GetEnable()) continue;
		double frac = GetFraction(i);
		sumFrac+=frac;
		outStr<<std::setw(10)<<resoList.at(i)->GetName()<<":    "<<frac<<"\n";
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

unsigned int AmpSumIntensity::GetNumberOfResonances() { return resoList.size(); }

int AmpSumIntensity::GetIdOfResonance(std::string name){
	for(unsigned int i=0; i< resoList.size(); i++)
		if(resoList.at(i)->GetName()==name) return i;
	return -999;
}

std::string AmpSumIntensity::GetNameOfResonance(unsigned int id){
	if(id < 0 || id > resoList.size() )
		throw std::runtime_error("AmpSumIntensity::getNameOfResonance() | "
				"Invalid resonance ID! Resonance not found?");
	return resoList.at(id)->GetName();
}

double AmpSumIntensity::GetMagnitude(std::string name) {
	int id = GetIdOfResonance(name);
	return GetMagnitude(id);
}

double AmpSumIntensity::GetMagnitude(unsigned int id) {
	if(id < 0 || id > resoList.size() )
		throw std::runtime_error("AmpSumIntensity::getAmpMagnitude() | "
				"Invalid resonance ID! Resonance not found?");
	return resoList.at(id)->GetMagnitude();
}

double AmpSumIntensity::GetPhase(std::string name) {
	int id = GetIdOfResonance(name);
	return GetPhase(id);
}

double AmpSumIntensity::GetPhase(unsigned int id) {
	if(id < 0 || id > resoList.size() )
		throw std::runtime_error("AmpSumIntensity::getAmpPhase() | "
				"Invalid resonance ID! Resonance not found?");
	return resoList.at(id)->GetPhase();
}

double AmpSumIntensity::GetIntegral(std::string name) {
	int id = GetIdOfResonance(name);
	return GetIntegral(id);
}

double AmpSumIntensity::GetIntegral(unsigned int id) {
	if(id < 0 || id > resoList.size() )
		throw std::runtime_error("AmpSumIntensity::getAmpIntegral() | "
				"Invalid resonance ID! Resonance not found?");
	//hard code normalization here to save cpu time
	//return 1;
	//recalculate integral for each call
	return resoList.at(id)->totalIntegral();
}

double AmpSumIntensity::GetFraction(std::string name) {
	int id = GetIdOfResonance(name);
	return GetFraction(id);
}

double AmpSumIntensity::GetFraction(unsigned int id) {
	if(id < 0 || id > resoList.size() )
		throw std::runtime_error("AmpSumIntensity::getAmpFraction() | "
				"Invalid resonance ID! Resonance not found?");
	return resoList.at(id)->GetMagnitude()/integral();
}
bool AmpSumIntensity::GetEnable(std::string name) {
	int id = GetIdOfResonance(name);
	return GetEnable(id);
}

bool AmpSumIntensity::GetEnable(unsigned int id) {
	if(id < 0 || id > resoList.size() )
		throw std::runtime_error("AmpSumIntensity::GetEnable() | "
				"Invalid resonance ID! Resonance not found?");
	return resoList.at(id)->GetEnable();
}

std::shared_ptr<AmpAbsDynamicalFunction>
AmpSumIntensity::GetResonance(std::string name) {
	int id = GetIdOfResonance(name);
	return GetResonance(id);
}

std::shared_ptr<AmpAbsDynamicalFunction>
AmpSumIntensity::GetResonance(unsigned int id) {
	if(id < 0 || id > resoList.size() )
		throw std::runtime_error("AmpSumIntensity::getResonance() | "
				"Invalid resonance ID! Resonance not found?");
	return resoList.at(id);
}

double AmpSumIntensity::averageWidth(){
	double avWidth = 0;
	double sum = 0;
	for(int i=0; i<resoList.size(); i++){
		avWidth += std::norm(resoList.at(i)->GetMagnitude())*resoList.at(i)->GetWidth();
		sum += std::norm(resoList.at(i)->GetMagnitude());
	}
	avWidth /= sum;
	return avWidth;
}
