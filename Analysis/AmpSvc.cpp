/*
 * AmpSvc.cpp
 *
 *  Created on: Aug 14, 2013
 *      Author: weidenka
 */
// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include  <cstdlib>

// Physics Interface header files go here
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Physics/AmplitudeSum/AmplitudeSetup.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

using namespace std;

int main(int argc, char **argv){
//	std::cout.setstate(std::ios::failbit) ;
	if(argc!=9&& argc!=6){
		cout<<"Usage: get fcn value at m23 m13 -> "<<argv[0]<<" <modelFile> <M> <m1> <m2> <m3> <m23> <m13> <rndVal>"<<endl;
		cout<<"Usage: get max value -> "<<argv[0]<<" <modelFile> <M> <m1> <m2> <m3>"<<endl;
		return 0;
	}
	std::string inputModel = argv[1];
	double M = atof(argv[2]);
	double m1 = atof(argv[3]);
	double m2 = atof(argv[4]);
	double m3 = atof(argv[5]);
	double m23;
	double m13;
	double rndVal;
	if(argc==9){
		m23=atof(argv[6]);
		m13=atof(argv[7]);
		rndVal=atof(argv[8]);
	}

	//load resonances
	AmplitudeSetup ini(inputModel);
//	cout << "loaded file " << ini.getFileName() << " with " << ini.getResonances().size() << " resonances!" << std::endl;
	AmpSumIntensity testBW(M, 0.0, m1, m2, m3, ini);
//	testBW.printAmps();
	ParameterList par; testBW.fillStartParVec(par);
	if(argc==6) {
//		std::cout.clear() ;
		cout<<"rEsUlT "<<testBW.getMaxVal()<<endl;
	}
	else{
		std::vector<double> x;
		x.push_back(m23);
		x.push_back(m13);
		double eval = testBW.intensity(x,par);
//		std::cout.clear() ;
		cout<<"rEsUlT ";
		if(eval>rndVal) cout<<"1"<<endl;
		else cout<<"0"<<endl;

	}
}
