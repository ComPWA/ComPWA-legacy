//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
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

using ComPWA::ControlParameter;
using ComPWA::DataReader::RootReader;
using ComPWA::Physics::BreitWigner::BreitWigner;
using ComPWA::ParameterList;
using ComPWA::DataReader::Data;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt" << std::endl;
  std::cout << std::endl;

  string file="test/2Part-4vecs.root";
  //RootReader myReader("test/2Part-4vecs.root");
  shared_ptr<RootReader> myReader(new RootReader(file, "data",true));
//  shared_ptr<BreitWigner> testBW(new BreitWigner(0.,5.));
  shared_ptr<ComPWA::Amplitude> testBW(new BreitWigner(0.,5.));

  shared_ptr<ControlParameter> testEsti = ComPWA::Estimator::MinLogLH::MinLogLH::createInstance(testBW, myReader, std::shared_ptr<Data>(), std::shared_ptr<Data>());
  //if(!testEsti) return 0;
  ParameterList minPar;
  testBW->FillParameterList(minPar);
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

  shared_ptr<ControlParameter> testEsti2 = ComPWA::Estimator::ChiOneD::ChiOneD::createInstance(testBW, myReader);
  minPar.SetParameterValue(0,1.5);
  minPar.SetParameterValue(1,0.3);
  result = testEsti2->controlParameter(minPar);
  cout << "1dim Fit optimal Chi: " << result << endl;

  return 0;
}
