/*
 * DPKinematics.cpp
 *
 *  Created on: Oct 16, 2013
 *      Author: weidenka
 */

#include "Physics/DPKinematics/DPKinematics.hpp"
#include "Physics/DPKinematics/DPpoint.hpp"

void DPKinematics::init(){
	m23_sq_min= ((m2+m3)*(m2+m3));
	m23_sq_max=((M-m1)*(M-m1));
	m13_sq_min=((m1+m3)*(m1+m3));
	m13_sq_max=((M-m2)*(M-m2));
	m12_sq_min=((m1+m2)*(m1+m2));
	m12_sq_max=((M-m3)*(M-m3));
	m23_min=((m2+m3)); m23_max=((M-m1));
	m13_min=((m1+m3)); m13_max=((M-m2));
	m12_min=((m1+m2)); m12_max=((M-m3));
};
DPKinematics::DPKinematics(double _M, double _Br, double _m1, double _m2, double _m3, std::string _name1, std::string _name2, std::string _name3):
				M(_M), Br(_Br), m1(_m1), m2(_m2), m3(_m3), name1(_name1), name2(_name2), name3(_name3)
{
	init();
};
DPKinematics::DPKinematics(const DPKinematics& other):
				M(other.M), Br(other.Br),
				m1(other.m1), m2(other.m2), m3(other.m3),
				name1(other.name1), name2(other.name2), name3(other.name3)
{
	init();
};
bool DPKinematics::isWithinDP(DPpoint point) const{
	return 0;
}
bool DPKinematics::isWithinDP(double m13, double m23) const{
	return 0;
}
