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

class AmpKinematics{
public:
	enum barrierType{BWPrime,BW,none};
	AmpKinematics(DoubleParameter& mR,int subSys, int spin, int m, int n, barrierType type,
			DoubleParameter& mesonRadius,DoubleParameter& motherRadius);
	AmpKinematics(double ma, double mb, double mc, double M, DoubleParameter& mR, int subSys, barrierType type,
			int spin, int m, int n, DoubleParameter& mesonR, DoubleParameter& motherR);
	AmpKinematics(const AmpKinematics& other);
	virtual ~AmpKinematics(){};

	double getECMS3_12(){return 1;};
	double getECMS1_12(){return 1;};
	double EtoP(double e){return 1;};
	virtual void setDecayMasses(double, double, double, double);
	void setBarrierType(int) { };
	void setBarrierRadi(double mesonRadius, double motherRadius) {
		_mesonRadius.SetValue(mesonRadius); _motherRadius.SetValue(motherRadius);
	};

	static double qValue(double, double, double);
	static double FormFactor(double z0, double z,unsigned int spin);
	static double BlattWeiss(double x, double mR, double ma, double mb, double spin, double mesonRadius);

	double q0()  const { return qValue(_mR.GetValue(), _ma, _mb ); };
//	double q(double, double, double)  const;
	double q(double x)  const { return qValue(x, _ma, _mb); };
	double q(double x,double ma, double mb)  const { return qValue(x, ma, mb); };
	double BLres2(double x) const { return BlattWeiss(x,_mR.GetValue(),_ma,_mb,_spin,_mesonRadius.GetValue()); }
	double BLmother2(double x) const { return BlattWeiss(x,_mR.GetValue(),_ma,_mb,_spin,_motherRadius.GetValue()); }//is the spin correct here, or do we need to spin of mother particle?
//	double BLmother2(double x) const;

	double getResMass() {return _mR.GetValue();};
//	double getSpin() {return _spin;}; //needs to be declared in AmpAbsDynamicalFunction
	double getM() {return _m;};
	double getN() {return _n;};
	double getMesonRadius() {return _mesonRadius.GetValue();};
	double getMotherRadius() {return _motherRadius.GetValue();};

protected:
	double _ma, _mb, _mc, _M;
	DoubleParameter _mR;
	unsigned int _subSys;
	barrierType _type;
	unsigned int _spin; int _m; int _n;
	DoubleParameter _mesonRadius, _motherRadius;

private:
};

#endif /* AMPKINEMATICS_HPP_ */
