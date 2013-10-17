/*
 * DPpoint2.cpp
 *
 *  Created on: Oct 16, 2013
 *      Author: weidenka
 */


#include "Physics/DPKinematics/DPKinematics.hpp"
#include "Physics/DPKinematics/DataPoint.hpp"

dataPoint* dataPoint::inst = NULL;

double dataPoint::getMsq(int subsys){
	switch(subsys){
	case 3:
//		m12.setVal(1);
		return m12;break;
	case 4:
		return m13;break;
	case 5:
		return m23;break;
	default:
		std::cout<<"DPpoint2 getMsq(): wrong subsys!"<<std::endl;
		//			return RooRealVar();
	}
	return -999;
};
void dataPoint::setMsq(int subsys, double val){
	switch(subsys){
	case 3:
//		m12.setVal(val);
		m12=val; break;
	case 4:
//		m13.setVal(val);
		m13=val;break;
	case 5:
//		m23.setVal(val);
		m23=val;break;
	default:
		std::cout<<"DPpoint2 setMsq(): wrong subsys!"<<std::endl;
	}
	return;
};
