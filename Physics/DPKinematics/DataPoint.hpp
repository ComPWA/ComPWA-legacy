/*
 * DPpoint2.hpp
 *
 *  Created on: Oct 16, 2013
 *      Author: weidenka
 */
#ifndef DPPOINT2_HPP_
#define DPPOINT2_HPP_
#include "RooRealVar.h"

class DPKinematics;

class dataPoint
{
private:
	dataPoint(DPKinematics kin) : DPKin(kin) {};
	dataPoint(dataPoint const&){};
	dataPoint(){};

	static dataPoint* inst;
public:
	static dataPoint* instance(DPKinematics kin){
		if(!inst) inst = new dataPoint(kin);
		inst->setKinematics(kin);
		return inst;
	};
	static dataPoint* instance(){
		if(!inst){
			std::cout<<"dataPoint: No instance created yet. Create instance using dataPoint::instance(DPKinematics kin) first!"<<std::endl;
			exit(1);
		}
		return inst;
	};
	bool isWithinDP() const{ return 0; };

	double getM(int subsys) {return sqrt(getMsq(subsys));};
	double getM(int a, int b) {return getM(a+b);};//daughter1 and daughter2
	double getMsq(int subsys);
	double getMsq(int a, int b) {return getM(a+b);};//daughter1 and daughter2
	//	RooRealVar& getMspec(int subsys);
	//	RooRealVar& getMspec(int a, int b) {return getMspec(a+b);};//daughter1 and daughter2
	void setKinematics(DPKinematics kin) { DPKin=kin; };
	void setM(int sys, double val) { setMsq(sys, val*val); };
	void setM(int a, int b, double val) { setM(a+b, val);};
	void setMsq(int, double);
	void setMsq(int a, int b, double val) { setMsq(a+b, val);};
	//	void setM(std::string, std::string, double);
	DPKinematics DPKin;

protected:
	//	RooRealVar m13;
	//	RooRealVar m23;
	//	RooRealVar m12;
	double m13;
	double m23;
	double m12;
};


#endif /* DPPOINT2_HPP_ */
