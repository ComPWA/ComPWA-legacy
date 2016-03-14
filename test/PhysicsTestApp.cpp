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

using ComPWA::Amplitude;
using ComPWA::Physics::BreitWigner::BreitWigner;
using ComPWA::ParameterList;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt" << std::endl;
  std::cout << std::endl;

  shared_ptr<Amplitude> testBW(new BreitWigner(0.,5.));
  vector<double> x;
  x.push_back(1.5);
  ParameterList par;
  //par.push_back(shared_ptr<PWAParameter>(new PWAGenericPar<double>(1.5, 0.5, 2.5, 0.1)));
 // par.push_back(shared_ptr<PWAParameter>(new PWAGenericPar<double>(0.3, 0.1, 0.5, 0.05)));
  testBW->copyParameterList(par);
  ParameterList intensL = testBW->intensity(x, par);
  double BWpdf = intensL.GetDoubleParameter(0)->GetValue();
  cout << "BreitWigner Intensity: " << BWpdf << endl;
  cout << "BreitWigner Integral: " << testBW->integral(par) << endl;
  cout << "BreitWigner Intensity normalized: " << BWpdf / testBW->integral(par) << endl;

  cout << "Done ..." << endl << endl;

  return 0;
}
