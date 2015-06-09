//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

#ifndef AMPKINEMATICS_HPP_
#define AMPKINEMATICS_HPP_
#include <cstdlib>
#include <math.h>
#include <iostream>
#include "Core/Parameter.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

class AmpKinematics{
public:
	enum barrierType{BWPrime,BW,none};
	AmpKinematics(std::shared_ptr<DoubleParameter> mR,int subSys, int spin, int m, int n,
			std::shared_ptr<DoubleParameter> mesonRadius, std::shared_ptr<DoubleParameter> motherRadius);
	virtual ~AmpKinematics(){};

	void setBarrierType(int) { };
	void setBarrierRadi(double mesonRadius, double motherRadius) {
		_mesonRadius->SetValue(mesonRadius); _motherRadius->SetValue(motherRadius);
	};

	/** Convert width of resonance to coupling
	 *
	 * Implementation of Eq.47-21 of PDG2014. Only valid for narrow, isolated resonances.
	 * @param mSq invariant mass
	 * @param mR mass of resonance
	 * @param width width of resonance in channel [a,b]
	 * @param ma mass of particle a
	 * @param mb mass of particle b
	 * @return
	 */
	static std::complex<double> widthToCoupling(double mSq, double mR, double width, double ma, double mb, double spin, double mesonRadius);
	/** Convert coupling to width
	 *
	 * Convert coupling to channel (@ma,@mb) to partial width. Only valid for narrow, isolated resonances.
	 * Implementation of inverted Eq.47-21 of PDG2014.
	 * @param mSq invariant mass
	 * @param mR mass of resonance
	 * Implementation of inverted Eq.47-21 of PDG2014. Only valid for narrow, isolated resonances.
	 * @param mSq invariant mass
	 * @param mR mass of resonance
	 * @param g coupling to channel [a,b]
	 * @param ma mass of particle a
	 * @param mb mass of particle b
	 * @return
	 */
	static std::complex<double> couplingToWidth(double mSq, double mR, double g, double ma, double mb, double spin, double mesonRadius);
	/** Calculate Break-up momentum squared
	 *
	 * Calculate Break-up momentum at energy @sqrtS for particles with masses @ma and @mb.
	 * From PDG2014 Eq.46-20a. Below threshold the function is analytically continued.
	 * @param sqrtS center-of-mass energy
	 * @param ma mass particle A
	 * @param mb mass particle B
	 * @return |break-up momentum|
	 */
	static double qSqValue(double sqrtS, double ma, double mb);
	/** Calculate Break-up momentum
	 *
	 * Calculate Break-up momentum at energy @sqrtS for particles with masses @ma and @mb.
	 * From PDG2014 Eq.46-20a. Below threshold the function is analytically continued.
	 * @param sqrtS center-of-mass energy
	 * @param ma mass particle A
	 * @param mb mass particle B
	 * @return |break-up momentum|
	 */
	static std::complex<double> qValue(double sqrtS, double ma, double mb);
	/** Two body phsp factor
	 *
	 * From PDG2014 Eqn.47-2
	 * @param sqrtS invariant mass of particles A and B
	 * @param ma Mass of particle A
	 * @param mb Mass of particle B
	 * @return
	 */
	static std::complex<double> phspFactor(double sqrtS, double ma, double mb);
	static double FormFactor(double sqrtS, double ma, double mb, double spin, double mesonRadius);

	double getResMass() {return _mR->GetValue();};
//	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	double getM() {return _m;};
	double getN() {return _n;};
	double getMesonRadius() {return _mesonRadius->GetValue();};
	double getMotherRadius() {return _motherRadius->GetValue();};

	inline virtual bool isSubSys(const unsigned int subSys) const{ return (subSys==_subSys); };
	virtual double evaluateWignerD(dataPoint& point) {
		if(_spin==0) return 1.0;//save some computing time
		return _wignerD.evaluate(point);
	};

protected:
	double _ma, _mb, _mc, _M;
	//! resonance mass
	std::shared_ptr<DoubleParameter> _mR;
	unsigned int _subSys;
	unsigned int _spin; int _m; int _n;
	std::shared_ptr<DoubleParameter> _mesonRadius, _motherRadius;
	AmpWigner2 _wignerD;

private:
};

#endif /* AMPKINEMATICS_HPP_ */
