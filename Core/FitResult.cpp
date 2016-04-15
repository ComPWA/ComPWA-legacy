/*
 * FitResult.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: weidenka
 */
#include "Core/FitResult.hpp"

void FitResult::writeText(std::string filename){
	std::ofstream myfile;
	myfile.open(filename);
	genOutput(myfile);
	myfile.close();
	return;
};

void FitResult::writeSimpleText(std::string filename){
	std::ofstream myfile;
	myfile.open(filename);
	genSimpleOutput(myfile);
	myfile.close();
	return;
};

double FitResult::shiftAngle(double v){
	double originalVal = v;
	double val = originalVal;
	double pi = PhysConst::instance()->getConstValue("Pi");
	while(val> pi) val-=2*pi;
	while(val< -pi ) val+=2*pi;
	if(val!=originalVal)
		BOOST_LOG_TRIVIAL(info) << "shiftAngle(): shifting parameter from "<<originalVal<< " to "<<val<<"!";
	return val;
}

void FitResult::genSimpleOutput(std::ostream& out)
{
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		std::shared_ptr<DoubleParameter> outPar = finalParameters.GetDoubleParameter(o);
		out<<outPar->GetValue()<<" "<<outPar->GetError()<<" ";
	}
	out<<"\n";

	return;
}

void FitResult::setFinalParameters(ParameterList finPars)
{
	finalParameters.DeepCopy(finPars);
}

void FitResult::printFitParameters(TableFormater* tableResult)
{
	bool printTrue=0, printInitial=0;
	if(trueParameters.GetNParameter()) printTrue=1;
	if(initialParameters.GetNParameter()) printInitial=1;
	unsigned int parErrorWidth = 22;
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++)
		if(finalParameters.GetDoubleParameter(o)->GetErrorType()==ErrorType::ASYM) parErrorWidth=33;

	tableResult->addColumn("Nr");
	tableResult->addColumn("Name",15);
	if(printInitial) tableResult->addColumn("Initial Value",parErrorWidth);
	tableResult->addColumn("Final Value",parErrorWidth);
	if(printTrue) tableResult->addColumn("True Value",10);
	if(printTrue) tableResult->addColumn("Deviation",9);
	tableResult->header();

	std::shared_ptr<DoubleParameter> iniPar, outPar, truePar;
	for(unsigned int o=0;o<finalParameters.GetNDouble();o++){
		try{
			outPar = finalParameters.GetDoubleParameter(o);
		} catch (BadParameter& bad){
			BOOST_LOG_TRIVIAL(error) << "FitResult::printFitParameters() | "
					"can't access parameter of final parameter list!";
			throw;
		}
		if(printInitial){
			try{
				iniPar = initialParameters.GetDoubleParameter(outPar->GetName());
			} catch (BadParameter& bad){
				BOOST_LOG_TRIVIAL(error) << "FitResult::printFitParameters() | "
						"can't access parameter '"<<outPar->GetName()<<"' of initial parameter list!";
				throw;
			}
		}
		if(printTrue){
			try{
				truePar = trueParameters.GetDoubleParameter(outPar->GetName());
			} catch (BadParameter& bad){
				BOOST_LOG_TRIVIAL(error) << "FitResult::printFitParameters() | "
						"can't access parameter '"<<outPar->GetName()<<"' of true parameter list!";
				*tableResult << "not found"<< " - ";
				throw;
			}
		}

		ErrorType errorType = outPar->GetErrorType();
		bool isFixed = outPar->IsFixed();
		bool isAngle=0, isMag=0;
		if(outPar->GetName().find("phase")!=string::npos) isAngle=1;//is our Parameter an angle?
		if(outPar->GetName().find("mag")!=string::npos) isMag=1;//is our Parameter an angle?
		if(isAngle && !isFixed) {
			outPar->SetValue( shiftAngle(outPar->GetValue()) ); //shift angle to the interval [-pi;pi]
			if(printInitial) iniPar->SetValue( shiftAngle(iniPar->GetValue()) );
			if(printTrue) truePar->SetValue( shiftAngle(truePar->GetValue()) );
		}
		if(isMag && !isFixed) {
			outPar->SetValue( std::fabs(outPar->GetValue()) ); //abs value of parameter is magnitude
			if(printInitial) iniPar->SetValue( std::fabs(iniPar->GetValue()) );
			if(printTrue) truePar->SetValue( std::fabs(truePar->GetValue()) );
		}

		*tableResult << o << outPar->GetName();
		if(printInitial) *tableResult << *iniPar;// |nr.| name| inital value|
		if(isFixed) *tableResult<<" ";
		else
			*tableResult << *outPar;//final value
		if(printTrue){
			*tableResult << *truePar;
			double pi = PhysConst::instance()->getConstValue("Pi");
			double pull = (truePar->GetValue()-outPar->GetValue() );
			if(isAngle && !isFixed) { //shift pull by 2*pi if that reduces the deviation
				while( pull<0 && pull<-pi) pull+=2*pi;
				while( pull>0 && pull>pi) pull-=2*pi;
			}
			if(outPar->HasError()){
				if( errorType == ErrorType::ASYM && pull < 0)
					pull /= outPar->GetErrorLow();
				else if( errorType == ErrorType::ASYM && pull > 0)
					pull /= outPar->GetErrorHigh();
				else
					pull /= outPar->GetError();
			}
			if( !std::isnan(pull) )
				*tableResult << pull;
			else
				*tableResult << " ";
		}
	}
	tableResult->footer();

	return;
}
void FitResult::printFitFractions(TableFormater* tab)
{
	BOOST_LOG_TRIVIAL(info) << " FitResult::printFitFractions() | "
			"Calculating fit fractions!";
	auto itrAmp = _ampVec.begin();
	for( ; itrAmp!=_ampVec.end(); ++itrAmp){
		printFitFractions(tab, (*itrAmp));
	}
}

void FitResult::printFitFractions(TableFormater* fracTable,
		std::shared_ptr<Amplitude> amp)
{
	ParameterList ffList;
	calcFraction(ffList, amp);
	double sum=0, sumErrorSq=0;

	fracTable->Reset();

	std::string ampName = amp->GetName();
	//print matrix
	fracTable->addColumn("Resonance: "+ampName,40);//add empty first column
	fracTable->addColumn("Fraction",15);//add empty first column
	fracTable->addColumn("Error",15);//add empty first column
	fracTable->header();
	for(unsigned int i=0;i<ffList.GetNDouble(); ++i){
		std::shared_ptr<DoubleParameter> tmpPar = ffList.GetDoubleParameter(i);
		std::string resName = tmpPar->GetName();

		//Remove amplitude name from string
		std::string::size_type strPos = resName.find(ampName);
		if (strPos != std::string::npos)
			resName.erase(strPos, ampName.length());

		*fracTable << resName
			<< tmpPar->GetValue()
			<< tmpPar->GetError(); //assume symmetric errors here
		sum += tmpPar->GetValue();
		sumErrorSq += tmpPar->GetError()*tmpPar->GetError();
	}
	fracTable->delim();
	*fracTable << "Total" << sum << sqrt(sumErrorSq) ;
	fracTable->footer();

	return;
}

void FitResult::calcFraction()
{
	if(!fractionList.GetNDouble()) {
		calcFraction(fractionList);
		calcFractionError();
	} else
		BOOST_LOG_TRIVIAL(warning) << "FitResult::calcFractions() | "
				"Fractions already calculated. Skip!";
}

resonanceItr findResonancePartner(std::shared_ptr<Amplitude> amp, resonanceItr res)
{
	auto name = (*res)->GetName();
	auto it = amp->GetResonanceItrFirst();
	for( ; it != amp->GetResonanceItrLast(); ++it){ //fill matrix
		if( it == res ) continue;
		auto name2 = (*it)->GetName();
		if(name2.find(name) != std::string::npos) return it;
	}
	return res;

}

void FitResult::calcFraction(ParameterList& parList, std::shared_ptr<Amplitude> amp)
{
	double norm = 1.0;
	std::string ampName = amp->GetName();

	/* Unbinned efficiency correction in the FunctionTree does not provide
	 * an integral w/o efficiency correction. We have to calculate it here.
	 */
	try{
		norm = amp->GetIntegral();
	} catch (std::exception& ex){
		BOOST_LOG_TRIVIAL(error)<< "FitResult::calcFraction() | "
				"Normalization can't be calculated: "<<ex.what();
		throw;
	}

	BOOST_LOG_TRIVIAL(debug)<<"FitResult::calcFraction() | "
			"Amplitude "<<ampName<< " Norm="<<norm;

	//Start loop over resonances
	auto it = amp->GetResonanceItrFirst();
	for( ; it != amp->GetResonanceItrLast(); ++it){ //fill matrix
		if( (*it)->GetName().find("_CP")!=std::string::npos ) continue;

		// We search for a partner resonance and add it to the integral
		auto it2 = findResonancePartner(amp, it);

		// GetIntegralInterference returns the integal Int( A*B+B*A ),
		// including the complex coefficienct
		double nom = amp->GetIntegralInterference(it,it);
		if( it != it2 ){// Int |A+B|^2 = |A|^2 + |B|^2 + A*B + B*A
			double tmp22 = amp->GetIntegralInterference(it2,it2);
			double tmp12 = amp->GetIntegralInterference(it,it2);
			BOOST_LOG_TRIVIAL(debug) << "FitResult::calcFraction() | Calculating"
					<<" amplitude integral for composed amplitudes "
					<<(*it)->GetName()<<" and "<<(*it2)->GetName()<<": "
					<<"(11) "<<nom <<" (22) "<<tmp22 <<" (12) "<<tmp12
					<<" Total: "<<nom+tmp22+tmp12;
			nom += tmp22;
			nom += tmp12;
		} else {
			BOOST_LOG_TRIVIAL(debug) << "FitResult::calcFraction() | Resonance "
					"integal for "<<(*it)->GetName()<<": "<<nom;
		}

		std::string resName = ampName+" "+(*it)->GetName()+"_FF";
		std::shared_ptr<DoubleParameter> magPar = (*it)->GetMagnitudePar();
		double mag = magPar->GetValue(); //value of magnitude
		double magError = 0;
		if(magPar->HasError())
			magError = magPar->GetError(); //error of magnitude

		parList.AddParameter(
				std::shared_ptr<DoubleParameter>(
						new DoubleParameter(
								resName,
//								mag*mag*resInt/norm,
								nom/norm,
								std::fabs(2*(nom/mag)/norm * magError)
						)
				)
		);
	}

}


void FitResult::calcFraction(ParameterList& parList)
{
	if( !_ampVec.size() )
		throw std::runtime_error("FitResult::calcFractions() | "
				"No amplitude set, can't calculate fractions!");

	if(parList.GetNDouble())
		throw std::runtime_error("FitResult::calcFractions() | "
				"ParameterList not empty!");

	//	_amp->UpdateParameters(finalParameters); //update parameters in amplitude
	double norm =-1;

	//Start loop over amplitudes
	auto ampItr = _ampVec.begin();
	for( ; ampItr != _ampVec.end(); ++ampItr){
		calcFraction(parList, (*ampItr));
	}

	return;
}
