//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//--------------------------------------------------------------------------------

#include "HelicityIntensity.hpp"

namespace HelicityFormalism {

HelicityIntensity::HelicityIntensity() {
	// TODO Auto-generated constructor stub

}

HelicityIntensity::~HelicityIntensity() {
	// TODO Auto-generated destructor stub
}


const double HelicityIntensity::integral() {

}
const double HelicityIntensity::integral(ParameterList& par) {

}
const double HelicityIntensity::normalization() {

}
const double HelicityIntensity::normalization(ParameterList& par) {

}
double HelicityIntensity::getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen) {

}
double HelicityIntensity::getMaxVal(std::shared_ptr<Generator> gen) {

}


const ParameterList& HelicityIntensity::intensity(std::vector<double> point,
		ParameterList& par) {
	setParameterList(par);
	dataPoint dataP(point);
	return intensity(dataP);
}

const ParameterList& HelicityIntensity::intensity(dataPoint& point,
		ParameterList& par) {
	setParameterList(par);
	return intensity(point);
}

const ParameterList& HelicityIntensity::intensityNoEff(dataPoint& point) {
	std::complex<double> intensity = 0;
	if (Kinematics::instance()->isWithinPhsp(point)) {
	  //convert the data point into the boosted 4 vectors!


		for (unsigned int i = 0; i < amplitude_trees_.size(); ++i) {
			intensity += amplitude_trees_[i].evaluate(kinematics_trees_[i]);
		}
		intensity = pow(std::abs(intensity), 2.0);
	}

	if (intensity != intensity) {
		BOOST_LOG_TRIVIAL(error)<<"Intensity is not a number!!";
		intensity = 0;
	}
	result.SetParameterValue(0, intensity);
	return result;
}

const ParameterList& HelicityIntensity::intensity(dataPoint& point) {
	intensityNoEff(point);
	double eff = efficiency_->evaluate(point);
	result.SetParameterValue(0, result.GetDoubleParameter(0)->GetValue() * eff);
	return result;
}

const bool HelicityIntensity::fillStartParVec(ParameterList& outPar){
	outPar = ParameterList(params_);
	return true;
}

void HelicityIntensity::setParameterList(ParameterList& par){
	//parameters varied by Minimization algorithm
	if(par.GetNDouble()!=params_.GetNDouble())
		throw std::runtime_error("setParameterList(): size of parameter lists don't match");
	for(unsigned int i=0; i<params_.GetNDouble(); i++){
		std::shared_ptr<DoubleParameter> p = params_.GetDoubleParameter(i);
		if(!p->IsFixed()){
			p->SetValue(par.GetDoubleParameter(i)->GetValue());
			p->SetError(par.GetDoubleParameter(i)->GetError());
		}
	}
	return;
}

void HelicityIntensity::printAmps() {

}
void HelicityIntensity::printFractions() {

}

Amplitude* HelicityIntensity::Clone() {

}

} /* namespace HelicityFormalism */
