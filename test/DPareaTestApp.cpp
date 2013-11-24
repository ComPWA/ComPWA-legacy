//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff
//-------------------------------------------------------------------------------
//! Test App for calculation of DP area via numeric monte-carlo integration
/*!
 * @file DPareaTestApp.cpp
 * The app calculates the phase space are of a particle X with a mass of 1GeV decaying to 3 gamma's
 * From geometrical considerations the expected phasespace is a symmetric triangle with edge length of 1.0GeV.
 * Therefore the expected phase space size is A=1/2 * 1.0 * 1.0!
 */

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>


// Physics Interface header files go here
#include "Physics/DPKinematics/DalitzKinematics.hpp"
#include "Physics/DPKinematics/DataPoint.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){

	DalitzKinematics* kin = DalitzKinematics::createInstance("J/psi","gamma","gamma","gamma");
	double area = kin->getDParea();

	std::cout<<"DPareaApp: Phase space area expected from the decay X-> gamma gamma gamma with M(X)=1.0GeV "\
			"is a symmetric triangle with edge length of 1.0GeV. The expected area is therefore A = 0.5*1.0*1.0"<<std::endl;
	std::cout<<"DPareaApp: calculated phase space area: "<<area<<std::endl;
	return 0;
}
