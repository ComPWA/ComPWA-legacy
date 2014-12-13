#include "Physics/Background/FlatBackground.hpp"
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "gsl/gsl_monte.h"
#include "gsl/gsl_monte_plain.h"
#include "gsl/gsl_monte_miser.h"
#include "gsl/gsl_monte_vegas.h"

const ParameterList& FlatBackground::intensity(dataPoint& point){
	double AMPpdf=0;
	if(Kinematics::instance()->isWithinPhsp(point)) AMPpdf=1;
	double eff=eff_->evaluate(point);
	result.SetParameterValue(0,AMPpdf*eff);
	return result;
}
const ParameterList& FlatBackground::intensity(std::vector<double> point, ParameterList& par){
	setParameterList(par);
	dataPoint dataP; dataP.setVal(0,point[0]); dataP.setVal(1,point[1]);
	return intensity(dataP);
}
const ParameterList& FlatBackground::intensityNoEff(dataPoint& point) {
	result.SetParameterValue(0,1.);
	return result;
}

std::shared_ptr<FunctionTree> FlatBackground::setupBasicTree(
		allMasses& theMasses,allMasses& toyPhspSample,std::string suffix)
{
	BOOST_LOG_TRIVIAL(debug) << "FlatBackground::setupBasicTree() generating new tree!";
	if(theMasses.nEvents==0){
		BOOST_LOG_TRIVIAL(error) << "FlatBackground::setupBasicTree() data sample empty!";
		return std::shared_ptr<FunctionTree>();
	}
	if(toyPhspSample.nEvents==0){
		BOOST_LOG_TRIVIAL(error) << "FlatBackground::setupBasicTree() toy sample empty!";
		return std::shared_ptr<FunctionTree>();
	}
//	double _dpArea = Kinematics::instance()->getPhspVolume();
	//------------Setup Tree---------------------
	std::shared_ptr<FunctionTree> newTree = std::shared_ptr<FunctionTree>(new FunctionTree());
//	//------------Setup Tree Pars---------------------
//	std::shared_ptr<MultiDouble> m23sq = std::shared_ptr<MultiDouble>( new MultiDouble("m23sq",theMasses.masses_sq.at( std::make_pair(2,3) )) );
//	std::shared_ptr<MultiDouble> m13sq = std::shared_ptr<MultiDouble>( new MultiDouble("m13sq",theMasses.masses_sq.at( std::make_pair(1,3) )) );
//	std::shared_ptr<MultiDouble> m12sq = std::shared_ptr<MultiDouble>( new MultiDouble("m12sq",theMasses.masses_sq.at( std::make_pair(1,2) )) );
//	std::shared_ptr<MultiDouble> m23sq_phsp = std::shared_ptr<MultiDouble>( new MultiDouble("m23sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(2,3) )) );
//	std::shared_ptr<MultiDouble> m13sq_phsp = std::shared_ptr<MultiDouble>( new MultiDouble("m13sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(1,3) )) );
//	std::shared_ptr<MultiDouble> m12sq_phsp = std::shared_ptr<MultiDouble>( new MultiDouble("m12sq_phsp",toyPhspSample.masses_sq.at( std::make_pair(1,2) )) );
//	//----Strategies needed
//	std::shared_ptr<MultAll> mmultStrat = std::shared_ptr<MultAll>(new MultAll(ParType::MCOMPLEX));
	std::shared_ptr<MultAll> mmultDStrat = std::shared_ptr<MultAll>(new MultAll(ParType::MDOUBLE));
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

	std::shared_ptr<MultiDouble> one =std::shared_ptr<MultiDouble>(
			new MultiDouble("one", std::vector<double>(theMasses.nEvents,1.)) );
	newTree->createHead("Amplitude"+suffix, mmultDStrat);
	newTree->createLeaf("oneVect",one, "Amplitude"+suffix);
	newTree->createLeaf("one",1., "Amplitude"+suffix);
	return newTree;
}

double FlatBackground::evaluateEff(double x[], size_t dim) {
	if(dim!=2) return 0;
	//set data point: we assume that x[0]=m13 and x[1]=m23
	dataPoint point; point.setVal(1,x[0]); point.setVal(0,x[1]);
	if( !Kinematics::instance()->isWithinPhsp(point) ) return 0;//only integrate over phase space
	ParameterList res = intensity(point);
	double intens = *res.GetDoubleParameter(0);
	return intens;
}
double evaluateWrapperFlatBackground(double* x, size_t dim, void* param) {
	/* We need a wrapper here because intensity() is a member function of AmpAbsDynamicalFunction
	 * and can therefore not be referenced. But gsl_monte_function expects a function reference.
	 * As third parameter we pass the reference to the current instance of AmpAbsDynamicalFunction
	 */
	return static_cast<FlatBackground*>(param)->evaluateEff(x,dim);
};

const double FlatBackground::normalization(){
	size_t dim=2;
	double res=0.0, err=0.0;

	//set limits: we assume that x[0]=m13sq and x[1]=m23sq
	DalitzKinematics* kin = dynamic_cast<DalitzKinematics*>(Kinematics::instance());
	double xLimit_low[2] = {kin->m13_sq_min,kin->m23_sq_min};
	double xLimit_high[2] = {kin->m13_sq_max,kin->m23_sq_max};
	gsl_rng_env_setup ();
	const gsl_rng_type *T = gsl_rng_default; //type of random generator
	gsl_rng *r = gsl_rng_alloc(T); //random generator
	gsl_monte_function G = {&evaluateWrapperFlatBackground,dim, const_cast<FlatBackground*> (this)};

	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
	gsl_monte_vegas_integrate (&G, xLimit_low, xLimit_high, 2, _nCalls, r,s,&res, &err);
	gsl_monte_vegas_free(s);
	BOOST_LOG_TRIVIAL(info)<<"FlatBackground::integrate() Integration result for background: "<<res<<"+-"<<err<<" relAcc [%]: "<<100*err/res;

	return res;
}
