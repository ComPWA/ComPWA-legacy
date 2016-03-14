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
	virtual bool copyParameterList(ParameterList& par) { };

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
	virtual double getIntValue(std::string var1, double min1, double max1,
			std::string var2, double min2, double max2) { };

	const bool fillStartParVec(ParameterList& outPar);
	void setParameterList(ParameterList& par);

	void printAmps();
	void printFractions();

	Amplitude* Clone();
};

} /* namespace HelicityFormalism */

#endif /* HELICITYINTENSITY_HPP_ */
