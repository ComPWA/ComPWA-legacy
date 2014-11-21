/*
 * FitResult.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: weidenka
 */

#include <numeric>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/random/normal_distribution.hpp>
using namespace boost::log;

#include "Optimizer/Minuit2/MinuitResult.hpp"

#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

int multivariateGaussian(const gsl_rng *rnd, const int vecSize, const gsl_vector *in, const gsl_matrix *cov, gsl_vector *res){
	gsl_matrix *tmpM= gsl_matrix_alloc(vecSize,vecSize);
	gsl_matrix_memcpy(tmpM,cov);
	gsl_linalg_cholesky_decomp(tmpM);
	for(unsigned int i=0; i<vecSize; i++)
		gsl_vector_set( res, i, gsl_ran_ugaussian(rnd) );

	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, tmpM, res);
	gsl_vector_add(res,in);
	gsl_matrix_free(tmpM);

	return 0;
}


void MinuitResult::init(FunctionMinimum min){
	nRes = 0;
	MnUserParameterState minState = min.UserState();
	MnUserCovariance minuitCovMatrix = minState.Covariance();
	//	std::vector<double> minuitCovM = minState.Covariance().Data();//Covariance matrix is empty !?
	using namespace boost::numeric::ublas;

	if(!minState.HasCovariance()){
		BOOST_LOG_TRIVIAL(error)<<"MinuitResult: no valid correlation matrix available!";
		return;
	}
	/* Size of Minuit covariance vector is given by dim*(dim+1)/2.
	 * dim is the dimension of the covariance matrix.
	 * The dimension can therefore be calculated as dim = -0.5+-0.5 sqrt(8*size+1)
	 */
	unsigned int dim = minuitCovMatrix.Nrow();
	globalCC = minState.GlobalCC().GlobalCC();
	symmetric_matrix<double,upper> covMatrix(dim,dim);
	symmetric_matrix<double,upper> corrMatrix(dim,dim);
	//	if(minuitCovM.size()==dim*(dim+1)/2){
	for (unsigned i = 0; i < covMatrix.size1 (); ++ i)
		for (unsigned j = i; j < covMatrix.size2 (); ++ j){
			double entry = minuitCovMatrix(j,i);
			covMatrix (i, j) = entry;
			if(i==j) variance.push_back(sqrt(entry));
		}
	for (unsigned i = 0; i < covMatrix.size1 (); ++ i)
		for (unsigned j = i; j < covMatrix.size2 (); ++ j){
			double denom = variance[i]*variance[j];
			corrMatrix(i,j) = covMatrix(i,j)/denom;
		}
	//	} else BOOST_LOG_TRIVIAL(error)<<"MinuitResult: no valid correlation matrix available!";
	cov=covMatrix;
	corr=corrMatrix;
	initialLH = -1;
	finalLH = minState.Fval();
	edm= minState.Edm();
	isValid = min.IsValid();
	covPosDef = min.HasPosDefCovar();
	hasValidParameters = min.HasValidParameters();
	hasValidCov = min.HasValidCovariance();
	hasAccCov = min.HasAccurateCovar();
	hasReachedCallLimit = min.HasReachedCallLimit();
	hesseFailed = min.HesseFailed();
	errorDef = min.Up();
	nFcn = min.NFcn();

	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	useCorrelatedErrors=0;
	return;

}

void MinuitResult::setAmplitude(std::shared_ptr<Amplitude> newAmp){
	_amp=newAmp;
	nRes=_amp->getNumberOfResonances();
	return;
}

void  MinuitResult::smearParameterList(ParameterList& newParList){
	unsigned int nFree = 0;
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++)
		if(!finalParameters.GetDoubleParameter(o)->IsFixed()) nFree++;
	unsigned int t=0;
	gsl_vector* oldPar = gsl_vector_alloc(nFree);
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> outPar = finalParameters.GetDoubleParameter(o);
		if(outPar->IsFixed()) continue;
		gsl_vector_set(oldPar,t,outPar->GetValue());
		t++;
	}
	gsl_matrix* gslCov = gsl_matrix_alloc(nFree,nFree);
	for(unsigned int i=0; i<cov.size1();i++)
		for(unsigned int j=0; j<cov.size2();j++)
			gsl_matrix_set(gslCov,i,j,cov(i,j));

	gsl_vector* newPar = gsl_vector_alloc(nFree);
	multivariateGaussian(r, nFree,oldPar,gslCov,newPar);//generate set of smeared parameters

	newParList = ParameterList(finalParameters); //deep copy of finalParameters
	t=0;
	for(unsigned int o=0;o<newParList.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> outPar = newParList.GetDoubleParameter(o);
		if(outPar->IsFixed()) continue;
		outPar->SetValue(newPar->data[t]);//set floating values to smeard values
		t++;
	}
}
void MinuitResult::calcFractionError(std::vector<double>& fracError){
	nRes=_amp->getNumberOfResonances();
	if(useCorrelatedErrors){
		/* Exact error calculation */
		unsigned int numberOfSets = 1000;
		BOOST_LOG_TRIVIAL(info) << "Calculating errors of fit fractions using "<<numberOfSets<<" sets of parameters...";
		std::vector<std::vector<double> > fracVect;
		for(unsigned int i=0; i<numberOfSets; i++){
			ParameterList newPar; smearParameterList(newPar);
			//			std::cout<<i <<" ------ "<<newPar<<std::endl;
			_amp->setParameterList(newPar);//smear all free parameters according to cov matrix
			std::vector<double> tmp;
			calcFraction(tmp);
			fracVect.push_back(tmp);
			//			std::cout<<"adsfasd ";
			//			for(unsigned int i=0; i<tmp.size();i++) std::cout<<tmp.at(i)<<" ";
			//			std::cout<<std::endl;
		}
		//Calculate standart deviation
		for(unsigned int o=0;o<nRes;o++){
			double mean=0, sqSum=0., stdev=0;
			for(unsigned int i=0; i<fracVect.size();i++){
				double tmp = fracVect.at(i).at(o);
				mean += tmp;
				sqSum += tmp*tmp;
			}
			unsigned int s = fracVect.size();
			sqSum /= s;
			mean /= s;
			stdev = std::sqrt(sqSum - mean*mean); //this is crosscecked with the RMS of the distribution
			fracError.push_back(stdev);
		}
		//		std::cout<<"frac error: ";
		//			for(unsigned int i=0; i<fracError.size();i++) std::cout<<fracError.at(i)<<" ";
		//			std::cout<<std::endl;
		_amp->setParameterList(finalParameters); //set correct fit result
	} else {
		/* As a simple case we assume that the normalization integral doesn't have an error */
		BOOST_LOG_TRIVIAL(info) << "Calculating errors of fit fractions assuming that parameters "
				"are uncorrelated and that neglecting the error from normalization!";
		double norm = 0.0;
		if(!_amp->hasTree()) norm = _amp->integral();
		else {//if we have a tree, use it. Much faster especially in case of correlated errors in calcFractionError()
			_amp->getPhspTree()->recalculate();
			double phspVolume = Kinematics::instance()->getPhspVolume();
			/*We need the intensity over the PHSP without efficiency correction. Therefore we
			 * access node 'Amplitude' and sum up its values.*/
			//			std::shared_ptr<TreeNode> amplitudeNode = _amp->getPhspTree()->head()->getChildren().at(0)->getChildren().at(0);
			std::shared_ptr<TreeNode> amplitudeNode = _amp->getPhspTree()->head()->getChildNode("Amplitude");
			if(!amplitudeNode){
				BOOST_LOG_TRIVIAL(error)<<"MinuitResult::calcFractionError() : Can't find node 'Amplitude' in tree!";
				throw BadParameter("Node not found!");
			}

			if(amplitudeNode->getName()!="Amplitude") {
				BOOST_LOG_TRIVIAL(error)<<"MinuitResult::calcFractionError() : we expect node 'Amplitude' at that position,"
						" but found node "<<amplitudeNode->getName()<<". Probably the structure of the tree has changed!";
				throw BadParameter("Node not found!");
			}
			std::shared_ptr<MultiComplex> normPar = std::dynamic_pointer_cast<MultiComplex>(amplitudeNode->getValue());//node 'Amplitude'
			unsigned int numPhspEvents = normPar->GetNValues();
			for(unsigned int i=0; i<numPhspEvents;i++)
				norm+=abs(normPar->GetValue(i))*abs(normPar->GetValue(i));

			norm = norm*phspVolume/numPhspEvents; //correct calculation of normalization
		}

		for(unsigned int i=0;i<nRes; i++){ //fill matrix
			double resonanceInt = _amp->getTotalIntegral(i); //fit fraction of amplitude
			std::string parName = "mag_"+_amp->getNameOfResonance(i); //name of magnitude parameter
			//		std::shared_ptr<DoubleParameter> magPar = trueParameters.GetDoubleParameter(parName);
			std::shared_ptr<DoubleParameter> magPar = finalParameters.GetDoubleParameter(parName);
			double mag = magPar->GetValue(); //value of magnitude
			double magError = 0.0;
			if(magPar->HasError()) magError = magPar->GetError()->GetError();
			if(magPar->IsFixed()) fracError.push_back(0.);
			else fracError.push_back(2*mag*resonanceInt/norm * magError); // sigma_fraction = 2*|A| intRes/totalInt * sigma_A
			//			std::cout<<2*mag*resonanceInt/norm * magError<<std::endl;
		}
	}
	return;
}

void MinuitResult::calcFraction(std::vector<double>& frac){
	double norm = 1.36286;
	ParameterList currentPar;
	_amp->fillStartParVec(currentPar);
	if(!_amp->hasTree()) norm = _amp->integral();
	else {//if we have a tree, use it. Much faster especially in case of correlated errors in calcFractionError()
		_amp->getPhspTree()->recalculate();
		double phspVolume = Kinematics::instance()->getPhspVolume();
		/*We need the intensity over the PHSP without efficiency correction. Therefore we
		 * access node 'Amplitude' and sum up its values.*/
//		std::shared_ptr<TreeNode> amplitudeNode = _amp->getPhspTree()->head()->getChildren().at(0)->getChildren().at(0);
//		std::shared_ptr<TreeNode> amplitudeNode = _amp->getPhspTree()->head()->getChildNode("Amplitude");
		std::shared_ptr<TreeNode> amplitudeNode(_amp->getPhspTree()->head()->getChildNode("Amplitude"));
		if(!amplitudeNode){
			BOOST_LOG_TRIVIAL(error)<<"MinuitResult::calcFraction() : Can't find node 'Amplitude' in tree!";
			throw BadParameter("Node not found!");
		}
		if(amplitudeNode->getName()!="Amplitude") {
			BOOST_LOG_TRIVIAL(error)<<"MinuitResult::calcFraction() : we expect node 'Amplitude' at that position,"
					" but found node "<<amplitudeNode->getName()<<". Probably the structure of the tree has changed!";
			throw BadParameter("Node not found!");
		}

		std::shared_ptr<MultiComplex> normPar = std::dynamic_pointer_cast<MultiComplex>(amplitudeNode->getValue());//node 'Amplitude'
		unsigned int numPhspEvents = normPar->GetNValues();
		for(unsigned int i=0; i<numPhspEvents;i++)
			norm+=abs(normPar->GetValue(i))*abs(normPar->GetValue(i));

		norm = norm*phspVolume/numPhspEvents; //correct calculation of normalization
		//std::cout<<"Amplitude normalization: "<<norm<<std::endl;
		//		std::cout<<norm<<" "<<phspVolume<<" "<<numPhspEvents<<std::endl;
	}
	nRes=_amp->getNumberOfResonances();
	for(unsigned int i=0;i<nRes; i++){ //fill matrix
		double resonanceInt = 1.0;
		resonanceInt = _amp->getTotalIntegral(i);//this is simply the factor 2J+1, because the resonance is already normalized
		std::string parName = "mag_"+_amp->getNameOfResonance(i); //name of magnitude parameter
		std::shared_ptr<DoubleParameter> magPar = currentPar.GetDoubleParameter(parName);
		double mag = magPar->GetValue(); //value of magnitude
		frac.push_back(mag*mag*resonanceInt/norm); // f= |A|^2 * intRes/totalInt
	}
	return;
}

void MinuitResult::fractions(std::ostream& out){
	if(!_amp) {
		//		BOOST_LOG_TRIVIAL(error) << "MinuitResult::calculateFractions() | no amplitude set, can't calculate fractions!";
		return;
	}
	std::vector<double> fractions;
	BOOST_LOG_TRIVIAL(info) << "Calculating fit fractions...";
	calcFraction(fractions);
	std::vector<double> fractionsError;
	calcFractionError(fractionsError);
	double sum, sumErrorSq;

	//print matrix
	TableFormater fracTable(&out);
	fracTable.addColumn("Resonance",15);//add empty first column
	fracTable.addColumn("Fraction",15);//add empty first column
	fracTable.addColumn("Error",15);//add empty first column
	out<<"FIT FRACTIONS:"<<std::endl;
	fracTable.header();
	for(unsigned int i=0;i<fractions.size(); ++i){
		fracTable << _amp->getNameOfResonance(i) << fractions.at(i) << fractionsError.at(i);
		sum+=fractions.at(i);
		sumErrorSq+=fractionsError.at(i)*fractionsError.at(i);
	}
	fracTable.delim();
	fracTable << "Total" << sum << sqrt(sumErrorSq);
	fracTable.footer();

	return;
}
void MinuitResult::genOutput(std::ostream& out, std::string opt){
	bool printTrue=0;
	bool printParam=1, printCorrMatrix=1, printCovMatrix=1;
	if(opt=="P") {//print only parameters
		printCorrMatrix=0; printCovMatrix=0;
	}
	if(trueParameters.GetNParameter()) printTrue=1;
	TableFormater tableCov(&out);
	tableCov.addColumn(" ",15);//add empty first column
	out<<std::endl;
	out<<"--------------FIT RESULT----------------"<<std::endl;
	if(!isValid) out<<"		*** FIT RESULT NOT VALID! ***"<<std::endl;
	out<<"Initial Likelihood: "<<initialLH<<std::endl;
	out<<"Final Likelihood: "<<finalLH<<std::endl;
	out<<"Estimated distance to minimumn: "<<edm<<std::endl;
	out<<"Error definition: "<<errorDef<<std::endl;
	out<<"Number of calls: "<<nFcn<<std::endl;
	if(hasReachedCallLimit) out<<"		*** LIMIT OF MAX CALLS REACHED! ***"<<std::endl;
	out<<"CPU Time : "<<time/60<<"min"<<std::endl;
	out<<std::endl;

	if(!hasValidParameters) out<<"		*** NO VALID SET OF PARAMETERS! ***"<<std::endl;
	if(printParam){
		unsigned int parErrorWidth = 22;
		for(unsigned int o=0;o<finalParameters.GetNDouble();o++)
			if(finalParameters.GetDoubleParameter(o)->GetErrorType()==ErrorType::ASYM) parErrorWidth=33;

		out<<"PARAMETERS:"<<std::endl;
		TableFormater tableResult(&out);
		tableResult.addColumn("Nr");
		tableResult.addColumn("Name",15);
		tableResult.addColumn("Initial Value",parErrorWidth);
		tableResult.addColumn("Final Value",parErrorWidth);
		if(printTrue) tableResult.addColumn("True Value",13);
		if(printTrue) tableResult.addColumn("Deviation",13);
		tableResult.header();

		for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
			std::shared_ptr<DoubleParameter> iniPar = initialParameters.GetDoubleParameter(o);
			std::shared_ptr<DoubleParameter> outPar = finalParameters.GetDoubleParameter(o);
			ErrorType errorType = outPar->GetErrorType();
			bool isFixed = iniPar->IsFixed();
			bool isAngle=0;
			if(iniPar->GetName().find("phase")!=string::npos) isAngle=1;//is our Parameter an angle?
			if(isAngle && !isFixed) {
				outPar->SetValue( shiftAngle(outPar->GetValue()) ); //shift angle to the interval [-pi;pi]
			}

			tableResult << o << iniPar->GetName() << *iniPar ;// |nr.| name| inital value|
			if(isFixed) tableResult<<"FIXED";
			else {
				tableResult << *outPar;//final value
				tableCov.addColumn(iniPar->GetName(),15);//add columns in covariance matrix
			}
			if(printTrue){
				std::shared_ptr<DoubleParameter> truePar = trueParameters.GetDoubleParameter(iniPar->GetName());
				if(!truePar) {
					tableResult << "not found"<< " - ";
					continue;
				}
				tableResult << *truePar;
				double pull = (truePar->GetValue()-outPar->GetValue() );
				if( errorType == ErrorType::ASYM && pull < 0)
					pull /= outPar->GetError()->GetErrorLow();
				else if( errorType == ErrorType::ASYM && pull > 0)
					pull /= outPar->GetError()->GetErrorHigh();
				else
					pull /= outPar->GetError()->GetError();
				tableResult << pull;
			}
		}
		tableResult.footer();
	}

	if(!hasValidCov) out<<"		*** COVARIANCE MATRIX NOT VALID! ***"<<std::endl;
	if(!hasAccCov) out<<"		*** COVARIANCE MATRIX NOT ACCURATE! ***"<<std::endl;
	if(!covPosDef) out<<"		*** COVARIANCE MATRIX NOT POSITIVE DEFINITE! ***"<<std::endl;
	if(hesseFailed) out<<"		*** HESSE FAILED! ***"<<std::endl;
	if(hasValidCov){
		unsigned int n=0;
		if(printCovMatrix){
			out<<"COVARIANCE MATRIX:"<<std::endl;
			tableCov.header();
			for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
				std::shared_ptr<DoubleParameter> ppp = initialParameters.GetDoubleParameter(o);
				std::shared_ptr<DoubleParameter> ppp2 = finalParameters.GetDoubleParameter(o);
				if(ppp->IsFixed()) continue;
				tableCov << ppp->GetName();
				for(unsigned int t=0;t<cov.size1();t++) {
					if(n>=cov.size2()) { tableCov<< " "; continue; }
					if(t>=n)tableCov << cov(n,t);
					else tableCov << "";
				}
				n++;
			}
			tableCov.footer();
		}
		if(printCorrMatrix){
			out<<"CORRELATION MATRIX:"<<std::endl;
			TableFormater tableCorr(&out);
			tableCorr.addColumn(" ",15);//add empty first column
			tableCorr.addColumn("GlobalCC",10);//global correlation coefficient
			for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
				std::shared_ptr<DoubleParameter> ppp = finalParameters.GetDoubleParameter(o);
				if(ppp->IsFixed()) continue;
				tableCorr.addColumn(ppp->GetName(),15);//add columns in correlation matrix
			}
			tableCorr.header();
			n=0;
			for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
				std::shared_ptr<DoubleParameter> ppp = initialParameters.GetDoubleParameter(o);
				std::shared_ptr<DoubleParameter> ppp2 = finalParameters.GetDoubleParameter(o);
				if(ppp->IsFixed()) continue;
				tableCorr << ppp->GetName();
				//				if(globalCC.size()>o)
				tableCorr << globalCC[n]; //TODO: check if emtpy (don't know how this happened, but it did :)
				for(unsigned int t=0;t<corr.size1();t++) {
					if(n>=corr.size2()) { tableCorr<< " "; continue; }
					if(t>=n)tableCorr << corr(n,t);
					else tableCorr << "";
				}
				n++;
			}
			tableCorr.footer();
		}
	}
	fractions(out); //calculate and print fractions if amplitude is set

	return;
}
