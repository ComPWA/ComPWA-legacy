/*
 * FlatBackground.hpp
 *
 *  Created on: Dec 10, 2014
 *      Author: weidenka
 */

#ifndef FLATBACKGROUND_HPP_
#define FLATBACKGROUND_HPP_

#include "Core/Amplitude.hpp"
#include "Core/Kinematics.hpp"

class FlatBackground : public Amplitude
{
public:
	FlatBackground(std::shared_ptr<Efficiency> eff, unsigned int nCalls) : eff_(eff), _nCalls(nCalls) {
		result.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("AmpSumResult")));
	}

	//! wrapper function for amplitude value including efficiency
	double evaluateEff(double x[], size_t dim);
	virtual const double integral() { return Kinematics::instance()->getPhspVolume(); }
	virtual const double integral(ParameterList& par) {	return integral(); }
	virtual const double normalization();
	virtual const double normalization(ParameterList& par) { return normalization(); }

	virtual double getMaxVal(ParameterList& par, std::shared_ptr<Generator> gen) { return 1; }
	virtual double getMaxVal(std::shared_ptr<Generator> gen) { return 1; }

	virtual const ParameterList& intensity(dataPoint& point, ParameterList& par) { return intensity(point); }
	virtual const ParameterList& intensity(dataPoint& point);
	virtual const ParameterList& intensityNoEff(dataPoint& point);
	virtual const ParameterList& intensity(std::vector<double> point, ParameterList& par);

	virtual bool copyParameterList(ParameterList& outPar) { return 0; };
	virtual void setParameterList(ParameterList& par) {};

	//! Check of tree is available
	virtual bool hasTree(){
		return 1;
	}
	//! Getter function for function tree
	virtual std::shared_ptr<FunctionTree> getAmpTree(allMasses& theMasses,allMasses& toyPhspSample, std::string suffix=""){
		return setupBasicTree(theMasses,toyPhspSample, suffix);
	}
	//! Clone function
	virtual FlatBackground* Clone(){
		return (new FlatBackground(*this));
	}

	/* useless functions for this class */
	virtual void printAmps() {};
	virtual void printFractions() {};
	virtual double getIntValue(std::string st, double d1, double d2, std::string st2, double d3, double d4) { return 0.0; };

protected:
	/**Setup Basic Tree
	 *
	 * @param theMasses data sample
	 * @param toyPhspSample sample of flat toy MC events for normalization of the resonances
	 * @param opt Which tree should be created? "data" data Tree, "norm" normalization tree
	 * with efficiency corrected toy phsp sample or "normAcc" normalization tree with sample
	 * of accepted flat phsp events
	 */
	std::shared_ptr<FunctionTree> setupBasicTree(allMasses& theMasses,allMasses& toyPhspSample, std::string suffix="");

	std::shared_ptr<Efficiency> eff_;

	unsigned int _nCalls; //! precision for numeric integration
};




#endif /* FLATBACKGROUND_HPP_ */
