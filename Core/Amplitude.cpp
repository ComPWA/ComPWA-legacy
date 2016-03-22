/*
 * Amplitude.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: weidenka
 */



#include "Core/Amplitude.hpp"


void Amplitude::UpdateParameters(ParameterList& par)
{
	for(unsigned int i=0; i<params.GetNDouble(); i++){
		std::shared_ptr<DoubleParameter> p, pNew;
		p = params.GetDoubleParameter(i);
		try{
			pNew = par.GetDoubleParameter(p->GetName());
		} catch (std::exception& ex){
			BOOST_LOG_TRIVIAL(error) << "AmpSumIntensity::setParameterList() | "
					" Can not find parameter! "<<ex.what();
			throw;
		}
		//Update parameter
		p->UpdateParameter( pNew );
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
