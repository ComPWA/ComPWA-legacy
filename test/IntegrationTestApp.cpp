//! Test-Application for full fit with simple modules.
/*!
 * @file EstimatorTestApp.cpp
 * This tiny application tests a simple fit procedure with a set of simple
 * modules. It uses a simle \f$\chi^{2}\f$ estimator EIFChiOneD, it reads data
 * via the root-reader module DIFRootReader and uses the intensity provided by
 * the simple 1D-Breit-Wigner physics module PIFBW. The optimization of the
 * parameters is done with the Minuit2 module OIFMinuit. As result the
 * optimized parameters are printed to the terminal.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

// ComPWA header files go here
#include "DIFRootReader.hpp"
#include "PIFBW.hpp"
#include "EIFChiOneD.hpp"
#include "OIFMinuit.hpp"

//Test header files go here
#include "PolyFit.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  string file="test/2Part-4vecs.root";
  shared_ptr<DIFRootReader> myReader(new DIFRootReader(file));
  shared_ptr<PIFBW> testBW(new PIFBW());
  shared_ptr<EIFChiOneD> testEsti(new EIFChiOneD(testBW, myReader));
  shared_ptr<OIFBase> opti(new OIFMinuit(testEsti));

  // Initiate parameters
  double val[2], min[2], max[2], err[2];
  val[0] = 1.5; max[0] = 2.5; min[0] = 0.5; err[0] = 0.5;
  val[1] = 0.3; max[1] = 0.5; min[1] = 0.1; err[1] = 0.1;

  double genResult = opti->exec(2, val,  min, max, err);

  std::cout << "Minimized final par :\t" << genResult << std::endl;
  std::cout << "final M:\t" << val[0] << " +- " << err[0] << std::endl;
  std::cout << "final T:\t" << val[1] << " +- " << err[1] << std::endl;

  return 0;
}
