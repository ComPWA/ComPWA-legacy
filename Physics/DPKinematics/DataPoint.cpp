//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

#include "Physics/DPKinematics/DalitzKinematics.hpp"
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
		m12=val;
		break;
	case 4:
		m13=val;
		break;
	case 5:
		m23=val;
		break;
	default:
		std::cout<<"DPpoint2 setMsq(): wrong subsys!"<<std::endl;
	}
	return;
};
