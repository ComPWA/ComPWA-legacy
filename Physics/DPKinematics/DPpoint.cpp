/*
 * DPpoint.cpp
 *
 *  Created on: Oct 16, 2013
 *      Author: weidenka
 */


#include "Physics/DPKinematics/DPKinematics.hpp"
#include "Physics/DPKinematics/DPpoint.hpp"

DPpoint::DPpoint(DPKinematics kin) : DPKin(kin),
m13("m13","inv. mass 1 3", DPKin.m13_min,DPKin.m13_min),
m23("m23","inv. mass 2 3", DPKin.m23_min,DPKin.m23_min),
m12("m23","inv. mass 2 3", DPKin.m12_min,DPKin.m12_min)
{
};

RooRealVar& DPpoint::getM(int subsys){
	switch(subsys){
	case 3:
		m12.setVal(1);
		return m12;
	case 4:
		return m13;
	case 5:
		return m23;
	default:
		std::cout<<"DPpoint: wrong subsys!"<<std::endl;
		//			return RooRealVar();
	}
};

RooRealVar& DPpoint::getMspec(int subsys){
	switch(subsys){
	case 4: //are the cases correct?
		m12.setVal(1);
		return m12;
	case 5:
		return m13;
	case 3:
		return m23;
	default:
		std::cout<<"DPpoint: wrong subsys!"<<std::endl;
		//			return 0;
	}
};
