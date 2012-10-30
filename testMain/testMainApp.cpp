//============================================================================
// Name        : testMain.cpp
// Author      : Mathias Michel
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <memory>

#include "OptimizerTest.h"
#include "DataTest.h"
#include "PhysicsTest.h"

using namespace std;

int main(int argc, char **argv) {

	//========= c++11 test ==========
	int test=7;
	auto neu=test;
	shared_ptr<int> test_ptr(new int(8));
	auto* auto_ptr(new int(9));
	cout << "say hello to c++11" << endl << endl;
	cout << "auto type: " << neu << endl;
	cout << "shared pointer:  " << *test_ptr << endl;
	cout << "auto pointer:  " << *auto_ptr << endl;
	cout << endl << endl;

	//======= Optimizer Test =========
	cout << "Optimizer Interface test started" << endl << endl;
    //OptimizerTest(argc, argv);

	//======== Physics Test ==========
	cout << "Physics Interface test started" << endl << endl;
    PhysicsTest(argc, argv);

	//========== Data Test ===========
	cout << "Data Interface test started" << endl << endl;
    DataTest(argc, argv);

	//========== Finished? ===========
	cout << "Now check results please" << endl << endl;

	return 0;
}
