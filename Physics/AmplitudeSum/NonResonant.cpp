/*
 * NonResonant.cpp
 *
 *  Created on: Jan 13, 2015
 *      Author: weidenka
 */

#include "Physics/AmplitudeSum/NonResonant.hpp"


NonResonant::NonResonant(std::string name) : AmpAbsDynamicalFunction(name.c_str())
{
}

std::complex<double> NonResonant::dynamicalFunction(){
	return std::complex<double>(1,0);
}



