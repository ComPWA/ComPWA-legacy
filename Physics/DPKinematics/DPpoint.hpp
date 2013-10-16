/*
 * DPpoint.hpp
 *
 *  Created on: Oct 16, 2013
 *      Author: weidenka
 */
#ifndef DPPOINT_HPP_
#define DPPOINT_HPP_
#include "RooRealVar.h"
class DPKinematics;

class DPpoint
{
public:
	DPpoint(){};
	DPpoint(DPKinematics kin);

	bool isWithinDP() const{};

	RooRealVar& getM(int subsys);
	RooRealVar& getM(int a, int b) {return getM(a+b);};//daughter1 and daughter2
	RooRealVar& getMspec(int subsys);
	RooRealVar& getMspec(int a, int b) {return getMspec(a+b);};//daughter1 and daughter2
	//	RooRealVar& getM(std::string, std::string);
	//	void setM(int, double);
	//	void setM(int a, int b, double val) { setM(a+b, val);};
	//	void setM(std::string, std::string, double);
	const DPKinematics DPKin;
protected:
	RooRealVar m13;
	RooRealVar m23;
	RooRealVar m12;
};


#endif /* DPPOINT_HPP_ */
