/*
 * FitResult.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: weidenka
 */

#include <numeric>
#include <cmath>
#include <cstdlib>

#include <boost/archive/xml_oarchive.hpp>

#include "Core/ProgressBar.hpp"
#include "Core/Logging.hpp"
#include "Optimizer/Minuit2/MinuitResult.hpp"

MinuitResult::MinuitResult(std::shared_ptr<ControlParameter> esti,
		FunctionMinimum result) :
			useCorrelatedErrors(0), calcInterference(0)
			correlatedErrors_numberOfSets(200)
{
	std::shared_ptr<Estimator> est = std::static_pointer_cast<Estimator>(esti);
	_ampVec = est->getAmplitudes();
	penalty = est->calcPenalty();
	nEvents = est->getNEvents();
	init(result);
}

void MinuitResult::setResult(std::shared_ptr<ControlParameter> esti,
		FunctionMinimum result)
{
	std::shared_ptr<Estimator> est = std::static_pointer_cast<Estimator>(esti);
	_ampVec = est->getAmplitudes();
	penalty = est->calcPenalty();
	nEvents = est->getNEvents();
	init(result);
}

void MinuitResult::init(FunctionMinimum min){
	nRes = 0;
	MnUserParameterState minState = min.UserState();

	if(minState.HasCovariance()){
		MnUserCovariance minuitCovMatrix = minState.Covariance();
		/* Size of Minuit covariance vector is given by dim*(dim+1)/2.
		 * dim is the dimension of the covariance matrix.
		 * The dimension can therefore be calculated as
		 * dim = -0.5+-0.5 sqrt(8*size+1)
		 */
		nFreeParameter = minuitCovMatrix.Nrow();
		globalCC = minState.GlobalCC().GlobalCC();
		cov = std::vector<std::vector<double>>(
				nFreeParameter,std::vector<double>(nFreeParameter));
		corr = std::vector<std::vector<double>>(
				nFreeParameter,std::vector<double>(nFreeParameter));
		for (unsigned i = 0; i < nFreeParameter; ++i)
			for (unsigned j = i; j < nFreeParameter; ++j){
				cov.at(i).at(j) = minuitCovMatrix(j,i);
				cov.at(j).at(i) = minuitCovMatrix(j,i);//fill lower half
			}
		for (unsigned i = 0; i < nFreeParameter; ++i)
			for (unsigned j = i; j < nFreeParameter; ++j){
				corr.at(i).at(j) =
						cov.at(i).at(j) / sqrt( cov.at(i).at(i) *
								cov.at(j).at(j) );
				corr.at(j).at(i) = corr.at(i).at(j);//fill lower half
			}

	} else
		BOOST_LOG_TRIVIAL(error)
		<< "MinuitResult: no valid correlation matrix available!";
	initialLH = -1;
	finalLH = minState.Fval();
	edm= minState.Edm();
	isValid = min.IsValid();
	covPosDef = min.HasPosDefCovar();
	hasValidParameters = min.HasValidParameters();
	hasValidCov = min.HasValidCovariance();
	hasAccCov = min.HasAccurateCovar();
	hasReachedCallLimit = min.HasReachedCallLimit();
	edmAboveMax = min.IsAboveMaxEdm();
	hesseFailed = min.HesseFailed();
	errorDef = min.Up();
	nFcn = min.NFcn();

	//	Amplitude::FillAmpParameterToList(_ampVec, finalParameters);

	if( Amplitude::AmpHasTree(_ampVec) ) setUseTree(1);
	return;

}

void MinuitResult::genSimpleOutput(std::ostream& out){
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> outPar =
				finalParameters.GetDoubleParameter(o);
		out<<outPar->GetValue()<<" "<<outPar->GetError()<<" ";
	}
	out<<"\n";

	return;
}

void MinuitResult::setUseCorrelatedErrors(bool s, int nSets) {
	useCorrelatedErrors = s;
	if(nSets <= 0)
		throw std::runtime_error(
				"MinuitResult::setUseCorrelatedErrors() |"
				" nSets <=0. That makes no sense!"
				);
	correlatedErrors_numberOfSets = nSets;
	return;
}

void MinuitResult::calcFractionError()
{
	if(useCorrelatedErrors){/* Exact error calculation */
		BOOST_LOG_TRIVIAL(info) << "Calculating errors of fit fractions using "
				<<correlatedErrors_numberOfSets<<" sets of parameters...";

		//Setting up random number generator
		const gsl_rng_type * T;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		gsl_rng* rnd = gsl_rng_alloc (T);

		//convert to GSL objects
		gsl_vector* gslFinalPar = gsl_parameterList2Vec(finalParameters);
		gsl_matrix* gslCov = gsl_vecVec2Matrix(cov);
		gsl_matrix_print(gslCov); //DEBUG

		std::vector<ParameterList> fracVect;
		progressBar bar(correlatedErrors_numberOfSets);
		stringstream outFraction;
		for(unsigned int i=0; i<correlatedErrors_numberOfSets; i++){
			bar.nextEvent();
			gsl_vector* gslNewPar = gsl_vector_alloc(nFreeParameter);
			//generate set of smeared parameters
			multivariateGaussian( rnd, nFreeParameter,
					gslFinalPar, gslCov, gslNewPar );
			gsl_vector_print(gslNewPar);
			//deep copy of finalParameters
			ParameterList newPar = ParameterList(finalParameters);
			std::size_t t=0;
			for(std::size_t o=0;o<newPar.GetNDouble();o++){
				std::shared_ptr<DoubleParameter> outPar =
						newPar.GetDoubleParameter(o);
				if(outPar->IsFixed()) continue;
				//set floating values to smeared values
				outPar->SetValue(gslNewPar->data[t]);
				t++;
			}
			//free vector
			gsl_vector_free(gslNewPar);
			//update amplitude with smeared parameters
			Amplitude::UpdateAmpParameterList(_ampVec, newPar);
			ParameterList tmp;
			calcFraction(tmp);
			fracVect.push_back(tmp);

			/******* DEBUGGING *******/
			//			if(i==0){
			//				for(int t=0; t<newPar.GetNDouble(); t++){
			//					if( newPar.GetDoubleParameter(t)->IsFixed()) continue;
			//					outFraction << newPar.GetDoubleParameter(t)->GetName()<<":";
			//				}
			//				for(int t=0; t<tmp.GetNDouble(); t++)
			//					outFraction << tmp.GetDoubleParameter(t)->GetName()<<":";
			//				outFraction << "norm" << std::endl;
			//			}
			//			for(int t=0; t<newPar.GetNDouble(); t++){
			//				if( newPar.GetDoubleParameter(t)->IsFixed()) continue;
			//				outFraction << newPar.GetDoubleParameter(t)->GetValue()<<" ";
			//			}
			//			for(int t=0; t<tmp.GetNDouble(); t++)
			//				outFraction << tmp.GetDoubleParameter(t)->GetValue()<<" ";
			//			double norm = _amp->GetIntegral();
			//			outFraction << norm;
			//			outFraction << std::endl;
			/******* DEBUGGING *******/
		}
		BOOST_LOG_TRIVIAL(info)<<" ------- "<<outFraction.str();

		//free objects
		gsl_vector_free(gslFinalPar);
		gsl_matrix_free(gslCov);
		gsl_rng_free(rnd);

		int nRes=fractionList.GetNDouble();
		//Calculate standard deviation
		for(unsigned int o=0;o<nRes;o++){
			double mean=0, sqSum=0., stdev=0;
			for(unsigned int i=0; i<fracVect.size();i++){
				double tmp = fracVect.at(i).GetDoubleParameter(o)->GetValue();
				mean += tmp;
				sqSum += tmp*tmp;
			}
			unsigned int s = fracVect.size();
			sqSum /= s;
			mean /= s;
			//this is cross-checked with the RMS of the distribution
			stdev = std::sqrt(sqSum - mean*mean);
			fractionList.GetDoubleParameter(o)->SetError(stdev);
		}

		//Set correct fit result
		Amplitude::UpdateAmpParameterList(_ampVec, finalParameters);
	}
	return;
}

void MinuitResult::genOutput(std::ostream& out, std::string opt){
	bool printTrue=0;
	bool printParam=1, printCorrMatrix=1, printCovMatrix=1;
	if(opt=="P") {//print only parameters
		printCorrMatrix=0; printCovMatrix=0;
	}
	if(trueParameters.GetNParameter()) printTrue=1;
	out<<std::endl;
	out<<"--------------MINUIT2 FIT RESULT----------------"<<std::endl;
	if(!isValid) out<<"		*** MINIMUM NOT VALID! ***"<<std::endl;
	out<<std::setprecision(10);
	out<<"Initial Likelihood: "<<initialLH<<std::endl;
	out<<"Final Likelihood: "<<finalLH<<std::endl;

	out<<"Estimated distance to minimumn: "<<edm<<std::endl;
	if(edmAboveMax) out<<"		*** EDM IS ABOVE MAXIMUM! ***"<<std::endl;
	out<<"Error definition: "<<errorDef<<std::endl;
	out<<"Number of calls: "<<nFcn<<std::endl;
	if(hasReachedCallLimit)
		out<<"		*** LIMIT OF MAX CALLS REACHED! ***"<<std::endl;
	out<<"CPU Time : "<<time/60<<"min"<<std::endl;
	out<<std::setprecision(5)<<std::endl;


	if(!hasValidParameters)
		out<<"		*** NO VALID SET OF PARAMETERS! ***"<<std::endl;
	if(printParam){
		out<<"PARAMETERS:"<<std::endl;
		TableFormater* tableResult = new TableFormater(&out);
		printFitParameters(tableResult);
	}

	if(!hasValidCov)
		out<<"		*** COVARIANCE MATRIX NOT VALID! ***"<<std::endl;
	if(!hasAccCov)
		out<<"		*** COVARIANCE MATRIX NOT ACCURATE! ***"<<std::endl;
	if(!covPosDef)
		out<<"		*** COVARIANCE MATRIX NOT POSITIVE DEFINITE! ***"<<std::endl;
	if(hesseFailed)
		out<<"		*** HESSE FAILED! ***"<<std::endl;
	if(hasValidCov){
		unsigned int n=0;
		if(printCovMatrix){
			out<<"COVARIANCE MATRIX:"<<std::endl;
			TableFormater* tableCov = new TableFormater(&out);
			printCovarianceMatrix(tableCov);
		}
		if(printCorrMatrix){
			out<<"CORRELATION MATRIX:"<<std::endl;
			TableFormater* tableCorr = new TableFormater(&out);
			printCorrelationMatrix(tableCorr);
		}
	}
	out<<"FIT FRACTIONS:"<<std::endl;
	//calculate and print fractions if amplitude is set
	TableFormater tab(&out);
	printFitFractions(&tab);

	out<<std::setprecision(10);
	out<<"Final penalty term: "<<penalty<<std::endl;
	out<<"FinalLH w/o penalty: "<<finalLH-penalty<<std::endl;
	out<<"FinalLH w/ penalty: "<<finalLH<<std::endl;
	/* The Akaike (AIC) and Bayesian (BIC) information criteria are described in
	 * Schwarz, Anals of Statistics 6 No.2: 461-464 (1978)
	 * and
	 * IEEE Transacrions on Automatic Control 19, No.6:716-723 (1974) */
	out<<"AIC: "<<calcAIC()-penalty<<std::endl;
	out<<"BIC: "<<calcBIC()-penalty<<std::endl;
	double r=0;
	for(int i=0; i<fractionList.GetNDouble(); i++){
		double val = std::fabs(fractionList.GetDoubleParameter(i)->GetValue());
		if(val > 0.001) r++;
	}
	out<<"Number of Resonances > 10^-3: "<<r<<std::endl;

	if(calcInterference){
		auto ampItr = _ampVec.begin();
		for( ; ampItr != _ampVec.end(); ++ampItr)
			createInterferenceTable(out,(*ampItr));
	}

	out<<std::setprecision(5);//reset cout precision
	return;
}

void MinuitResult::createInterferenceTable(std::ostream& out,
		std::shared_ptr<Amplitude> amp)
{
		out<<"INTERFERENCE terms for "<<amp->GetName()<<": "<<std::endl;
		TableFormater* tableInterf = new TableFormater(&out);
		tableInterf->addColumn("Name 1",15);
		tableInterf->addColumn("Name 2",15);
		tableInterf->addColumn("Value",15);
		tableInterf->header();
		double sumInfTerms = 0;
		auto it = amp->GetResonanceItrFirst();
		for( ; it != amp->GetResonanceItrLast(); ++it){
			auto it2 = it;
			for( ; it2 != amp->GetResonanceItrLast(); ++it2){
				*tableInterf << (*it)->GetName();
				*tableInterf << (*it2)->GetName();
				double inf = amp->GetIntegralInterference(it,it2);
				*tableInterf << inf;
				sumInfTerms+=inf;
			}
		}
		tableInterf->delim();
		*tableInterf<<" "<<"Sum: "<<sumInfTerms;
		tableInterf->footer();
		out<<std::endl;
}

double MinuitResult::calcAIC()
{
	if(!fractionList.GetNDouble()) calcFraction();
	double r=0;
	for(int i=0; i<fractionList.GetNDouble(); i++){
		double val = fractionList.GetDoubleParameter(i)->GetValue();
		if(val > 0.001) r++;
	}
	return (finalLH+2*r);
}

double MinuitResult::calcBIC(){
	if(!fractionList.GetNDouble()) calcFraction();
	double r=0;
	for(int i=0; i<fractionList.GetNDouble(); i++){
		double val = fractionList.GetDoubleParameter(i)->GetValue();
		if(val > 0.001) r++;
	}
	return (finalLH+r*std::log(nEvents));
}

void MinuitResult::printCorrelationMatrix(TableFormater* tableCorr){
	if(!hasValidCov) return;
	tableCorr->addColumn(" ",15);//add empty first column
	tableCorr->addColumn("GlobalCC",10);//global correlation coefficient

	//add columns in correlation matrix
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> ppp =
				finalParameters.GetDoubleParameter(o);
		if(ppp->IsFixed()) continue;
		tableCorr->addColumn(ppp->GetName(),15);
	}

	unsigned int n=0;
	tableCorr->header();
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> ppp =
				finalParameters.GetDoubleParameter(o);
		if(ppp->IsFixed()) continue;
		*tableCorr << ppp->GetName();
		*tableCorr << globalCC.at(n);
		for(unsigned int t=0;t<corr.size();t++) {
			if(n>=corr.at(0).size()) { *tableCorr<< " "; continue; }
			if(t>=n)*tableCorr << corr.at(n).at(t);
			else *tableCorr << "";
		}
		n++;
	}
	tableCorr->footer();
	return;
}

void MinuitResult::printCovarianceMatrix(TableFormater* tableCov){
	if(!hasValidCov) return;
	tableCov->addColumn(" ",17);//add empty first column
	//add columns first
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		if(!finalParameters.GetDoubleParameter(o)->IsFixed())
			tableCov->addColumn(
					finalParameters.GetDoubleParameter(o)->GetName(),17);
	}

	unsigned int n=0;
	tableCov->header();
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> ppp =
				finalParameters.GetDoubleParameter(o);
		if(ppp->IsFixed()) continue;
		*tableCov << ppp->GetName();
		for(unsigned int t=0;t<cov.size();t++) {
			if(n>=cov.at(0).size()) { *tableCov<< " "; continue; }
			if(t>=n) *tableCov << cov.at(n).at(t);
			else *tableCov << "";
		}
		n++;
	}
	tableCov->footer();
	return;
}

void MinuitResult::writeXML(std::string filename){
	std::ofstream ofs(filename);
	boost::archive::xml_oarchive oa(ofs);
	oa << boost::serialization::make_nvp("FitParameters", finalParameters);
	oa << boost::serialization::make_nvp("FitFractions", fractionList);
	ofs.close();
	return;
}

void MinuitResult::writeTeX(std::string filename){
	std::ofstream out(filename);
	bool printTrue=0;
	if(trueParameters.GetNParameter()) printTrue=1;
	TableFormater* tableResult = new TexTableFormater(&out);
	printFitParameters(tableResult);
	if(hasValidCov){
		unsigned int n=0;
		TableFormater* tableCov = new TexTableFormater(&out);
		printCovarianceMatrix(tableCov);
		TableFormater* tableCorr = new TexTableFormater(&out);
		printCorrelationMatrix(tableCorr);
	}
	TableFormater* fracTable = new TexTableFormater(&out);
	//calculate and print fractions if amplitude is set
	printFitFractions(fracTable);
	out.close();
	return;
}

bool MinuitResult::hasFailed(){
	bool failed=0;
	if(!isValid) failed=1;
	//	if(!covPosDef) failed=1;
	//	if(!hasValidParameters) failed=1;
	//	if(!hasValidCov) failed=1;
	//	if(!hasAccCov) failed=1;
	//	if(hasReachedCallLimit) failed=1;
	//	if(hesseFailed) failed=1;

	return failed;
}
