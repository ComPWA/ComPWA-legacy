#include "Optimizer/Geneva/GenevaResult.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Geneva {

void GenevaResult::init(boost::shared_ptr<Gem::Geneva::GStartIndividual> min){
	finalLH = min->getBestCase();
	min->getPar(finalParameters);
    //ToDO: extract more info
	return;
}

void GenevaResult::genOutput(std::ostream& out, std::string opt){
	bool printTrue=0;
	if(trueParameters.GetNParameter()) printTrue=1;
	bool printParam=1;

	TableFormater tableCov(&out);
	tableCov.addColumn(" ",15);//add empty first column
	out<<std::endl;
	out<<"--------------GENEVA FIT RESULT----------------"<<std::endl;
	out<<"Final Likelihood: "<<finalLH<<std::endl;
	if(printParam){
		out<<"PARAMETERS:"<<std::endl;
		TableFormater* tableResult = new TableFormater(&out);
		printFitParameters(tableResult);
	}

	out<<"FIT FRACTIONS:"<<std::endl;
	TableFormater* fracTable = new TableFormater(&out);
	printFitFractions(fracTable); //calculate and print fractions if amplitude is set
	return;
}

} /* namespace Geneva */
} /* namespace Optimizer */
} /* namespace ComPWA */
