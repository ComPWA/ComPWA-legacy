//! Test-Application of the Physics-IF.
/*!
 * @file PhysicsTestApp.cpp
 * This tiny application tests the interface to the Physics-Module. The simple
 * implementation using a 1-dim Breit-Wigner is used and the intensity at the mean
 * of the distribution is printed.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

// Physics Interface header files go here
#include "Physics/BreitWigner/BreitWigner.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){

  shared_ptr<Amplitude> testBW(new BreitWigner(0.,5.));
  vector<double> x;
  x.push_back(1.5);
  ParameterList par;
  //par.push_back(shared_ptr<PWAParameter>(new PWAGenericPar<double>(1.5, 0.5, 2.5, 0.1)));
 // par.push_back(shared_ptr<PWAParameter>(new PWAGenericPar<double>(0.3, 0.1, 0.5, 0.05)));
  testBW->fillStartParVec(par);
  cout << "BreitWigner Intensity: " << testBW->intensity(x, par) << endl;
  cout << "BreitWigner Integral: " << testBW->integral(par) << endl;
  cout << "BreitWigner Intensity normalized: " << testBW->intensity(x, par) / testBW->integral(par) << endl;

  cout << "Done ..." << endl << endl;

  return 0;
}
