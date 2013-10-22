/*
 * AmpKinematics.hpp
 *
 *  Created on: Oct 12, 2013
 *      Author: weidenka
 */

#ifndef AMPKINEMATICS_HPP_
#define AMPKINEMATICS_HPP_
#include <cstdlib>
#include <math.h>
#include <iostream>
#include "Core/Parameter.hpp"

class AmpKinematics{
public:
	enum barrierType{BWPrime,BW,none};
	AmpKinematics(DoubleParameter,int, int, int, int, barrierType, double, double);
	AmpKinematics(double, double, double, double, DoubleParameter, int, barrierType, int, int, int, double, double);
	AmpKinematics(const AmpKinematics& other);
	~AmpKinematics(){};

	double getECMS3_12(){return 1;};
	double getECMS1_12(){return 1;};
	double EtoP(double e){return 1;};
	virtual void setDecayMasses(double, double, double, double);
	void setBarrierType(int) { };
	void setBarrierRadi(double mesonRadius, double motherRadius) {_mesonRadius=mesonRadius; _motherRadius=motherRadius; };

	virtual double integral() const = 0;
	double q0(double, double) const;
	double q0()  const { return q0( _ma, _mb ); };
	double q(double, double, double)  const;
	double q(double x)  const { return q(x, _ma, _mb); };
	double BLres2(double x) const;
	double BLmother2(double x) const;

	double getResMass() {return _mR.GetValue();};
//	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	double getM() {return _m;};
	double getN() {return _n;};
	double getMesonRadius() {return _mesonRadius;};
	double getMotherRadius() {return _motherRadius;};

protected:
	double _ma, _mb, _mc, _M;
	DoubleParameter _mR;
	barrierType _type;
	int _spin; int _m; int _n;
	int _subSys;
	double _mesonRadius, _motherRadius;

private:
	double FormFactor(double z0, double z) const;
};

#endif /* AMPKINEMATICS_HPP_ */
