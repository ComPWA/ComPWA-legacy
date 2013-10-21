#ifndef DPKINEMATICS_HPP_
#define DPKINEMATICS_HPP_

#include "RooAbsReal.h"

class DPpoint;
class dataPoint;

class DPKinematics
{
public:
	DPKinematics(){};
	void init();
	DPKinematics(double _M, double _Br, double _m1, double _m2, double _m3, std::string _name1, std::string _name2, std::string _name3);
	DPKinematics(const DPKinematics& other);
	double getThirdVariable(double, double) const;
	bool isWithinDP() const;
	bool isWithinDP(double m23, double m13, double m12=0) const;

	double lambda(double x, double y, double z)const;
	Double_t s2min(Double_t s1, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const;
	Double_t s2max(Double_t s1, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const;
	Double_t s3min(Double_t s1, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const;
	Double_t s3max(Double_t s1, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const;
	Double_t s1min(Double_t s2, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const;
	Double_t s1max(Double_t s2, Double_t m0, Double_t m1, Double_t m2, Double_t m3)const;
	Double_t s2min(Double_t s1)const { return s2min(s1,M,m1,m2,m3); };
	Double_t s2max(Double_t s1)const { return s2max(s1,M,m1,m2,m3); };
	Double_t s3min(Double_t s1)const { return s3min(s1,M,m1,m2,m3); };
	Double_t s3max(Double_t s1)const { return s3max(s1,M,m1,m2,m3); };
	Double_t s1min(Double_t s2)const { return s1min(s2,M,m1,m2,m3); };
	Double_t s1max(Double_t s2)const { return s1max(s2,M,m1,m2,m3); };

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
