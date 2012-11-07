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
#include "PIFBW.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){

  PIFBW* testBW = new PIFBW();
  cout << "BreitWigner Intensity: " << testBW->intensity(1.5, 1.5, 0.3) << endl;

  cout << "Done ..." << endl << endl;

  return 0;
}
