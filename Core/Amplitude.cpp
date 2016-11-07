/*
 * Amplitude.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: weidenka
 */

#include "Core/Amplitude.hpp"

namespace ComPWA {


void Amplitude::UpdateParameters(ParameterList& par)
{
	std::shared_ptr<DoubleParameter> pOld, pNew;

	/* First we check if new parameter list contains at least one matching
	 * parameter. Otherwise we skip! */
	int commonPar = 0;
	for(unsigned int i=0; i<params.GetNDouble(); i++){
		try{
			pNew = par.GetDoubleParameter(
					params.GetDoubleParameter(i)->GetName()
			);
		} catch (std::exception& ex){
			continue;
		}
		commonPar++;
	}
	if(commonPar == 0) return;

	/* If we have at least one matching parameter, we require that all
	 * parameters are contained in the new list */
	for(unsigned int i=0; i<params.GetNDouble(); i++){
		pOld = params.GetDoubleParameter(i);
		try{
			pNew = par.GetDoubleParameter( pOld->GetName() );
		} catch (std::exception& ex){
			BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::setParameterList() | "
					" Can not find parameter! "<<ex.what();
			throw;
		}
		//Update parameter
		pOld->UpdateParameter( pNew );
	}

	return;
}

void Amplitude::FillParameterList(ParameterList& outPar) const
{
	//Parameters are only added if they do not exist yet
	outPar.Append(params);
	outPar.RemoveDuplicates();
}

void Amplitude::FillAmpParameterToList(
		std::vector<std::shared_ptr<Amplitude> > ampVec,
		ParameterList& list)
{
	auto it = ampVec.begin();
	for( ; it!=ampVec.end(); ++it )
		(*it)->FillParameterList(list);

	list.RemoveDuplicates(); //Parameter should contain each enry only once
}

double Amplitude::intensityInterference(dataPoint& point,
		resonanceItr A, resonanceItr B)
{
	double intens = (
			(*A)->Evaluate(point)*std::conj((*B)->Evaluate(point))
	).real();
	if( A != B ) intens = 2*intens;

	return intens;
}

void Amplitude::SetAmpEfficiency(
		std::vector<std::shared_ptr<Amplitude> > ampVec,
		std::shared_ptr<Efficiency> eff)
{
	auto it = ampVec.begin();
	for( ; it!=ampVec.end(); ++it ){
		(*it)->SetEfficiency(eff);
	}
}

void Amplitude::UpdateAmpParameterList(
		std::vector<std::shared_ptr<Amplitude> > ampVec,
		ParameterList& list)
{
	auto it = ampVec.begin();
	for( ; it!=ampVec.end(); ++it )
		(*it)->UpdateParameters(list);
}

bool Amplitude::AmpHasTree(std::vector<std::shared_ptr<Amplitude> > ampVec)
{
	auto it = ampVec.begin();
	for( ; it!=ampVec.end(); ++it )
		if( !(*it)->hasTree() ) return 0;

	return 1;
}

void Amplitude::GetAmpFitFractions(
		std::vector<std::shared_ptr<Amplitude> > ampVec,
		ParameterList& parList)
{
	if( !ampVec.size() )
		throw std::runtime_error("FitResult::calcFractions() | "
				"No amplitude set, can't calculate fractions!");

	if(parList.GetNDouble())
		throw std::runtime_error("FitResult::calcFractions() | "
				"ParameterList not empty!");

	//	_amp->UpdateParameters(finalParameters); //update parameters in amplitude
	double norm =-1;

	//Start loop over amplitudes
	auto ampItr = ampVec.begin();
	for( ; ampItr != ampVec.end(); ++ampItr){
		(*ampItr)->GetFitFractions(parList);
	}

	return;
}

}
