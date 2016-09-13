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

HelicityIntensity::HelicityIntensity()
{
	// TODO Auto-generated constructor stub

}

HelicityIntensity::~HelicityIntensity()
{
	// TODO Auto-generated destructor stub
}


const double HelicityIntensity::GetIntegral()
{
	return 1.0;
}

const double HelicityIntensity::GetNormalization()
{
	return GetIntegral();
}

double HelicityIntensity::GetMaxVal(std::shared_ptr<Generator> gen)
{
	return 1.0;
}

const ParameterList& HelicityIntensity::intensity(dataPoint& point)
{
	intensityNoEff(point);
	double eff = efficiency_->evaluate(point);
	result.SetParameterValue(0, result.GetDoubleParameter(0)->GetValue() * eff);
	return result;
}

const ParameterList& HelicityIntensity::intensity(std::vector<double> point)
{
	dataPoint dataP(point);
	return intensity(dataP);
}

const ParameterList& HelicityIntensity::intensityNoEff(dataPoint& point)
{
	std::complex<double> intensity = 0;
	if (Kinematics::instance()->IsWithinPhsp(point)) {
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

void HelicityIntensity::to_str()
{

}


} /* namespace HelicityFormalism */
