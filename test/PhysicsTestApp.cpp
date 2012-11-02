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
