//! Test-Application to check c++11 support.
/*!
 * @file CX0TestApp.cpp
 * This tiny application test some features of the new c++11 standard. One can
 * compile and run this to check the c++11 support of the used compiler/system.
*/

#include <iostream>
#include <memory>

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

	return 0;
}
