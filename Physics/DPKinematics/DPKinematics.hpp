
/*
 * DPKinematics.hpp
 *
 *  Created on: Oct 18, 2013
 *      Author: weidenka
 *
 *		DPKinematics stores information about the dalitz plot kinematics, namly the
 *		masses. It can caluclate boundaries and check if a certain point lies
 *		within the kinematically allowed region.
 */
#ifndef DPKINEMATICS_HPP_
#define DPKINEMATICS_HPP_

#include <iostream>
class DPpoint;
class dataPoint;

class DPKinematics
{
public:
	DPKinematics(){};
	void init();
	DPKinematics(std::string _nameMother, std::string _name1, std::string _name2, std::string _name3);
	DPKinematics(double _M, double _Br, double _m1, double _m2, double _m3, std::string _name1, std::string _name2, std::string _name3);
	DPKinematics(const DPKinematics& other);
	double getThirdVariable(double, double) const;
	bool isWithinDP() const;
	bool isWithinDP(double m23, double m13, double m12=0) const;

	double lambda(double x, double y, double z)const;
	double s2min(double s1, double m0, double m1, double m2, double m3)const;
	double s2max(double s1, double m0, double m1, double m2, double m3)const;
	double s3min(double s1, double m0, double m1, double m2, double m3)const;
	double s3max(double s1, double m0, double m1, double m2, double m3)const;
	double s1min(double s2, double m0, double m1, double m2, double m3)const;
	double s1max(double s2, double m0, double m1, double m2, double m3)const;
	double s2min(double s1)const { return s2min(s1,M,m1,m2,m3); };
	double s2max(double s1)const { return s2max(s1,M,m1,m2,m3); };
	double s3min(double s1)const { return s3min(s1,M,m1,m2,m3); };
	double s3max(double s1)const { return s3max(s1,M,m1,m2,m3); };
	double s1min(double s2)const { return s1min(s2,M,m1,m2,m3); };
	double s1max(double s2)const { return s1max(s2,M,m1,m2,m3); };

	double M;
	double Br;
	double m1, m2, m3;
	std::string nameMother, name1, name2, name3;

	double m23_sq_min, m23_sq_max;
	double m13_sq_min, m13_sq_max;
	double m12_sq_min, m12_sq_max;
	double m23_min, m23_max;
	double m13_min, m13_max;
	double m12_min, m12_max;

};

#endif
