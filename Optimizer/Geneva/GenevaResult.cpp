#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
using namespace boost::log;

#include "Optimizer/Geneva/GenevaResult.hpp"


void GenevaResult::init(boost::shared_ptr<Gem::Geneva::GStartIndividual> min){
	finalLH = min->getBestKnownFitness();
    //ToDO: extract more info
	return;

}

void GenevaResult::genOutput(std::ostream& out){
	bool printTrue=0;
	if(trueParameters.GetNParameter()) printTrue=1;
	TableFormater tableCov(&out);
	tableCov.addColumn(" ",15);//add empty first column
	out<<std::endl;
	out<<"--------------FIT RESULT----------------"<<std::endl;
	//if(!isValid) out<<"		*** FIT RESULT NOT VALID! ***"<<std::endl;
	//out<<"Initial Likelihood: "<<initialLH<<std::endl;
	out<<"Final Likelihood: "<<finalLH<<std::endl;
	//out<<"Estimated distance to minimumn: "<<edm<<std::endl;
	//out<<"Error definition: "<<errorDef<<std::endl;
	//out<<"Number of calls: "<<nFcn<<std::endl;
	//if(hasReachedCallLimit) out<<"		*** LIMIT OF MAX CALLS REACHED! ***"<<std::endl;
	//out<<"CPU Time : "<<time/60<<"min"<<std::endl;
	//out<<std::endl;

	//if(!hasValidParameters) out<<"		*** NO VALID SET OF PARAMETERS! ***"<<std::endl;
	//out<<"PARAMETERS:"<<std::endl;
	//TableFormater tableResult(&out);
	//tableResult.addColumn("Nr");
	//tableResult.addColumn("Name",15);
	//tableResult.addColumn("Initial Value",20);
	//tableResult.addColumn("Final Value",20);
	//if(printTrue) tableResult.addColumn("True Value",10);
	//if(printTrue) tableResult.addColumn("true-final/error",16);
	//tableResult.header();

	/*for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> iniPar = initialParameters.GetDoubleParameter(o);
		std::shared_ptr<DoubleParameter> outPar = finalParameters.GetDoubleParameter(o);
		bool isFixed = iniPar->IsFixed();
		bool isAngle=0;
		if(iniPar->GetName().find("phi")!=string::npos) isAngle=1;//is our Parameter an angle?
		if(isAngle && !isFixed) outPar->SetValue( shiftAngle(outPar->GetValue()) ); //shift angle to the interval [-pi;pi]

		tableResult << o << iniPar->GetName() << *iniPar ;// |nr.| name| inital| value|
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
			tableResult << (truePar->GetValue()-outPar->GetValue() )/outPar->GetError();
		}
	}
	tableResult.footer();*/

	//if(!hasValidCov) out<<"		*** COVARIANCE MATRIX NOT VALID! ***"<<std::endl;
	//if(!hasAccCov) out<<"		*** COVARIANCE MATRIX NOT ACCURATE! ***"<<std::endl;
	//if(!covPosDef) out<<"		*** COVARIANCE MATRIX NOT POSITIVE DEFINITE! ***"<<std::endl;
	//if(!hesseFailed) out<<"		*** HESSE FAILED! ***"<<std::endl;
	/*if(hasValidCov){
		out<<"COVARIANCE MATRIX:"<<std::endl;
		tableCov.header();
		unsigned int n=0;
		for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
			std::shared_ptr<DoubleParameter> ppp = initialParameters.GetDoubleParameter(o);
			std::shared_ptr<DoubleParameter> ppp2 = finalParameters.GetDoubleParameter(o);
			if(ppp->IsFixed()) continue;
			tableCov << ppp->GetName();
			for(unsigned int t=0;t<cov.size1();t++) {
				if(n>=cov.size2()) { tableCov<< " "; continue; }
				tableCov << cov(n,t);
			}
			n++;
		}
		tableCov.footer();

		out<<"CORRELATION MATRIX:"<<std::endl;
		tableCov.header();
		n=0;
		for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
			std::shared_ptr<DoubleParameter> ppp = initialParameters.GetDoubleParameter(o);
			std::shared_ptr<DoubleParameter> ppp2 = finalParameters.GetDoubleParameter(o);
			if(ppp->IsFixed()) continue;
			tableCov << ppp->GetName();
			for(unsigned int t=0;t<corr.size1();t++) {
				if(n>=corr.size2()) { tableCov<< " "; continue; }
				tableCov << corr(n,t);
			}
			n++;
		}
		tableCov.footer();
	}*/

}
