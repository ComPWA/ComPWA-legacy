/*
 * FitResult.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: weidenka
 */



#include "Core/FitResult.hpp"
void FitResult::writeText(std::string filename){
	std::ofstream myfile;
	myfile.open(filename);
	genOutput(myfile);
	myfile.close();
	return;
};
void FitResult::writeSimpleText(std::string filename){
	std::ofstream myfile;
	myfile.open(filename);
	genSimpleOutput(myfile);
	myfile.close();
	return;
};
double FitResult::shiftAngle(double v){
	double val = v;
	double twoPi = 2*PhysConst::instance()->getConstValue("Pi");
	while(val> twoPi) val-=twoPi;
	while(val< -twoPi ) val+=twoPi;
	return val;
}
