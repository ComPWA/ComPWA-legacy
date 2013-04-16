//! Test-Application of a simple \f$\chi^{2}\f$ estimator.
/*!
 * @file EstimatorTestApp.cpp
 * This tiny application tests a simple \f$\chi^{2}\f$ estimator module. It reads data
 * via the root-reader module RootReader,hpp and uses the intensity provided by
 * the simple 1D-Breit-Wigner physics module BreitWigner.hpp. As result it prints a
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
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/BreitWigner/BreitWigner.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Estimator/ChiOneD/ChiOneD.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  string file="test/2Part-4vecs.root";
  //RootReader myReader("test/2Part-4vecs.root");
  shared_ptr<RootReader> myReader(new RootReader(file, true,"data"));
  shared_ptr<BreitWigner> testBW(new BreitWigner(0.,5.));

  shared_ptr<MinLogLH> testEsti(new MinLogLH(testBW, myReader));
  ParameterList minPar;
  testBW->fillStartParVec(minPar);
  double result=0;
  result = testEsti->controlParameter(minPar);
  cout << "1dim Fit optimal Likelihood: " << result << endl;

  minPar.SetParameterValue(0,1.7);
  minPar.SetParameterValue(1,0.3);
  result = testEsti->controlParameter(minPar);
  cout << "1.7 0.3: " << result << endl;

  minPar.SetParameterValue(0,1.3);
  minPar.SetParameterValue(1,0.3);
  result = testEsti->controlParameter(minPar);
  cout << "1.3 0.3: " << result << endl;

  minPar.SetParameterValue(0,1.5);
  minPar.SetParameterValue(1,0.2);
  result = testEsti->controlParameter(minPar);
  cout << "1.5 0.2: " << result << endl;

  minPar.SetParameterValue(0,1.5);
  minPar.SetParameterValue(1,0.4);
  result = testEsti->controlParameter(minPar);
  cout << "1.5 0.4:  " << result << endl;

  minPar.SetParameterValue(0,1.5);
  minPar.SetParameterValue(1,0.01);
  result = testEsti->controlParameter(minPar);
  cout << "1.5 0.01:  " << result << endl;

  shared_ptr<ChiOneD> testEsti2(new ChiOneD(testBW, myReader));
  minPar.SetParameterValue(0,1.5);
  minPar.SetParameterValue(1,0.3);
  result = testEsti2->controlParameter(minPar);
  cout << "1dim Fit optimal Chi: " << result << endl;

  return 0;
}
