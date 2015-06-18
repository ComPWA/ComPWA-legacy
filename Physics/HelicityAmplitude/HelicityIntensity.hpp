//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef HELICITYINTENSITY_HPP_
#define HELICITYINTENSITY_HPP_

#include "HelicityAmplitudeTree.hpp"

#include "Core/Amplitude.hpp"
#include "Core/Efficiency.hpp"

namespace HelicityFormalism {

class HelicityIntensity: public Amplitude {
	std::vector<HelicityAmplitudeTree> amplitude_trees_;
	std::vector<HelicityKinematicBoostTree> kinematics_trees_;
	std::shared_ptr<Efficiency> efficiency_;

	ParameterList params_;

public:
	HelicityIntensity();
	virtual ~HelicityIntensity();

	const double integral();
	const double integral(ParameterList& par);
	const double normalization();
	const double normalization(ParameterList& par);
	double getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen);
	double getMaxVal(std::shared_ptr<Generator> gen);

	const ParameterList& intensity(dataPoint& point, ParameterList& par);
	const ParameterList& intensity(dataPoint& point);
	const ParameterList& intensityNoEff(dataPoint& point);
	const ParameterList& intensity(std::vector<double> point,
			ParameterList& par);

	const bool fillStartParVec(ParameterList& outPar);
	void setParameterList(ParameterList& par);

	void printAmps();
	void printFractions();

	Amplitude* Clone();
};

} /* namespace HelicityFormalism */

#endif /* HELICITYINTENSITY_HPP_ */
