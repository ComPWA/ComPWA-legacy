//! Test-Application of a simple \f$\chi^{2}\f$ estimator.
/*!
 * @file EstimatorTestApp.cpp
 * This tiny application tests a simple \f$\chi^{2}\f$ estimator module. It reads data
 * via the root-reader module DIFRootReader,hpp and uses the intensity provided by
 * the simple 1D-Breit-Wigner physics module PIFBW.hpp. As result it prints a
 * calculated \f$\chi^{2}\f$ to the terminal.
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
#include "PWAParticle.hpp"
#include "PWAParameter.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  string file="test/2Part-4vecs.root";
  //DIFRootReader myReader("test/2Part-4vecs.root");
  shared_ptr<DIFRootReader> myReader(new DIFRootReader(file));
  shared_ptr<PIFBW> testBW(new PIFBW());

  shared_ptr<EIFChiOneD> testEsti(new EIFChiOneD(testBW, myReader));
  vector<PWAParameter<double> > minPar;
  minPar.push_back(PWAParameter<double>(1.5,1.,2.,0.1));
  minPar.push_back(PWAParameter<double>(0.3,0.1,0.2,0.01));
  double result=0;
  result = testEsti->controlParameter(minPar);
  cout << "1dim Fit optimal par Chi2: " << result << endl;

  cout << "Done ..." << endl << endl;
  return 0;
}
