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

	//! Clone function
	HelicityIntensity* Clone(std::string newName="") const {
		auto tmp = (new HelicityIntensity(*this));
		tmp->SetName(newName);
		return tmp;
	}

	const double GetIntegral();
	const double GetNormalization();
	double GetMaxVal(std::shared_ptr<Generator> gen);

	const ParameterList& intensity(dataPoint& point);
	const ParameterList& intensity(std::vector<double> point);
	const ParameterList& intensityNoEff(dataPoint& point);

	/** Integral value of amplitude in certain boundary
	 * Used for plotting a projection of a function in \p var1 in
	 * bin [\p min1, \p min2]. In this case we have to integrate over an
	 * other variable \p var2
	 * @param var1 first variable
	 * @param min1 minimal value of first variable
	 * @param max1 maximal value of first variable
	 * @param var2 second variable
	 * @param min2 minimal value of second variable
	 * @param max2 maximal value of second variable
	 * @return
	 */
	virtual double GetIntValue(std::string var1, double min1, double max1,
			std::string var2, double min2, double max2) { };

	void to_str();

};

} /* namespace HelicityFormalism */

#endif /* HELICITYINTENSITY_HPP_ */
