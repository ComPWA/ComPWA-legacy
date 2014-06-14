/*
 * FitResult.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: weidenka
 */


#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
using namespace boost::log;

#include "Optimizer/Minuit2/MinuitResult.hpp"


void MinuitResult::init(FunctionMinimum min){
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
	return;

}

void MinuitResult::genSimpleOutput(std::ostream& out){
	//	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
	//		std::shared_ptr<DoubleParameter> outPar = finalParameters.GetDoubleParameter(o);
	//		out<<outPar->GetName()<<"/D:"<<outPar->GetName()<<"err/D";
	//		if(o!=finalParameters.GetNDouble()-1) out<<":";
	//	}
	//	out<<"\n";
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> outPar = finalParameters.GetDoubleParameter(o);
		out<<outPar->GetValue()<<" "<<outPar->GetError()->GetError()<<" ";
	}
	out<<"\n";

	return;
}
void MinuitResult::setAmplitude(std::shared_ptr<Amplitude> newAmp){ _amp=newAmp; }

void MinuitResult::fractions(std::ostream& out){
	if(!_amp) {
		//		BOOST_LOG_TRIVIAL(error) << "MinuitResult::calculateFractions() | no amplitude set, can't calculate fractions!";
		return;
	}
	unsigned int nRes = _amp->getNumberOfResonances();
	std::cout<<nRes<<std::endl;
	fracError = boost::numeric::ublas::matrix<double>(nRes,2);
	double norm = _amp->integral();
	double sum = 0;
	double sumErrorSq = 0;
	unsigned int pos = 0;
	for(unsigned int i=0;i<nRes; i++){ //fill matrix
		double resonanceInt = _amp->getTotalIntegral(i); //fit fraction of amplitude
		std::string parName = "mag_"+_amp->getNameOfResonance(i); //name of magnitude parameter
		std::shared_ptr<DoubleParameter> magPar = finalParameters.GetDoubleParameter(parName);
		double mag = magPar->GetValue(); //value of magnitude
		double magError;
		if(magPar->IsFixed()) magError = 0.0;
		else{
			magError = variance.at(pos*2); //error on magnitude
			pos++;//we have to skip fixed variables
		}
		fracError(i,0) = mag*mag*resonanceInt/norm; // f= |A|^2 * intRes/totalInt
		fracError(i,1) = 2*mag*resonanceInt/norm*magError; // sigma_fraction = 2*|A| intRes/totalInt * sigma_A
		sum += fracError(i,0);
		sumErrorSq += fracError(i,1)*fracError(i,1);
		std::cout<<parName<<" "<<mag<<"+-"<<magError<<" "<<resonanceInt<<" "<<norm<<std::endl;
	}
	//print matrix
	TableFormater fracTable(&out);
	fracTable.addColumn("Resonance",15);//add empty first column
	fracTable.addColumn("Fraction",15);//add empty first column
	fracTable.addColumn("Error",15);//add empty first column
	out<<"FIT FRACTIONS:"<<std::endl;
	fracTable.header();
	for(unsigned int i=0;i<nRes; ++i)
		fracTable << _amp->getNameOfResonance(i) << fracError(i,0) << fracError(i,1);
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
		out<<"PARAMETERS:"<<std::endl;
		TableFormater tableResult(&out);
		tableResult.addColumn("Nr");
		tableResult.addColumn("Name",15);
		tableResult.addColumn("Initial Value",20);
		tableResult.addColumn("Final Value",30);
		if(printTrue) tableResult.addColumn("True Value",10);
		if(printTrue) tableResult.addColumn("true-final/error",16);
		tableResult.header();

		for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
			std::shared_ptr<DoubleParameter> iniPar = initialParameters.GetDoubleParameter(o);
			std::shared_ptr<DoubleParameter> outPar = finalParameters.GetDoubleParameter(o);
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
				tableResult << (truePar->GetValue()-outPar->GetValue() )/ *outPar->GetError();
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
				if(globalCC.size()>o)
					tableCorr << globalCC[o]; //TODO: check if emtpy (don't know how this happened, but it did :)
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
