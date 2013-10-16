#ifndef DPKINEMATICS_HPP_
#define DPKINEMATICS_HPP_

#include "RooAbsReal.h"

class DPpoint;

class DPKinematics
{
public:
	DPKinematics(){};
	void init();
	DPKinematics(double _M, double _Br, double _m1, double _m2, double _m3, std::string _name1, std::string _name2, std::string _name3);
	DPKinematics(const DPKinematics& other);
	bool isWithinDP(DPpoint point) const;
	bool isWithinDP(double m13, double m23) const;
	double M;
	double Br;
	double m1, m2, m3;
	std::string name1, name2, name3;

	double m23_sq_min, m23_sq_max;
	double m13_sq_min, m13_sq_max;
	double m12_sq_min, m12_sq_max;
	double m23_min, m23_max;
	double m13_min, m13_max;
	double m12_min, m12_max;

};

#endif
