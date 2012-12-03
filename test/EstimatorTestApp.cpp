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
#include "EIFMinLogLH.hpp"
#include "EIFChiOneD.hpp"
#include "PWAParticle.hpp"
#include "PWAParameter.hpp"
#include "PWAGenericPar.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  string file="test/2Part-4vecs.root";
  //DIFRootReader myReader("test/2Part-4vecs.root");
  shared_ptr<DIFRootReader> myReader(new DIFRootReader(file, true));
  shared_ptr<PIFBW> testBW(new PIFBW(0.,5.));

  shared_ptr<EIFMinLogLH> testEsti(new EIFMinLogLH(testBW, myReader));
  vector<shared_ptr<PWAParameter> > minPar;
  testBW->fillStartParVec(minPar);
  double result=0;
  result = testEsti->controlParameter(minPar);
  cout << "1dim Fit optimal Likelihood: " << result << endl;

  minPar[0]->SetValue(1.7);
  minPar[1]->SetValue(0.3);
  result = testEsti->controlParameter(minPar);
  cout << "1.7 0.3: " << result << endl;

  minPar[0]->SetValue(1.3);
  minPar[1]->SetValue(0.3);
  result = testEsti->controlParameter(minPar);
  cout << "1.3 0.3: " << result << endl;

  minPar[0]->SetValue(1.5);
  minPar[1]->SetValue(0.2);
  result = testEsti->controlParameter(minPar);
  cout << "1.5 0.2: " << result << endl;

  minPar[0]->SetValue(1.5);
  minPar[1]->SetValue(0.4);
  result = testEsti->controlParameter(minPar);
  cout << "1.5 0.4:  " << result << endl;

  minPar[0]->SetValue(1.5);
  minPar[1]->SetValue(0.01);
  result = testEsti->controlParameter(minPar);
  cout << "1.5 0.01:  " << result << endl;

  shared_ptr<EIFChiOneD> testEsti2(new EIFChiOneD(testBW, myReader));
  minPar[0]->SetValue(1.5);
  minPar[1]->SetValue(0.3);
  result = testEsti2->controlParameter(minPar);
  cout << "1dim Fit optimal Chi: " << result << endl;

  return 0;
}
