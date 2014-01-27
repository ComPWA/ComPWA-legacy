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
//! Test-Application for fit with simple modules using RunManager.
/*!
 * @file RunManagerTestApp.cpp
 * This tiny application tests a simple fit procedure with a set of simple
 * modules. It uses a simle \f$\chi^{2}\f$ estimator ChiOneD, it reads data
 * via the root-reader module RootReader and uses the intensity provided by
 * the simple 1D-Breit-Wigner physics module BreitWigner. The optimization of the
 * parameters is done with the Minuit2 module MinuitIF. The class RunManager
 * is used to perform the actual fitting procedure.
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
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/RunManager.hpp"

#include "Physics/DPKinematics/RootEfficiency.hpp"
//Test header files go here
#include "PolyFit.hpp"

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  boost::log::core::get()->set_filter(trivial::severity >= trivial::debug); //setting log level
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt" << std::endl;
  std::cout << std::endl;

  std::string file="test/2Part-4vecs.root";
  std::cout << "Load Modules" << std::endl;
  std::shared_ptr<Data> myReader(new RootReader(file, false,"data",true));
  std::shared_ptr<Amplitude> testBW(new BreitWigner(0.,5.));
  // Initiate parameters
  ParameterList par;
  par.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BWPos",1.7,0.5,2.5,0.1)));
  par.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("BWWidth",0.2,0.1,0.2,0.01)));
  std::shared_ptr<ControlParameter> testEsti = MinLogLH::createInstance(testBW, myReader); //TODO: <- should be done by runManager
  std::shared_ptr<Optimizer> opti(new MinuitIF(testEsti,par));
//  std::shared_ptr<RunManager> run(new RunManager(myReader, testEsti, testBW, opti, eff));
  std::shared_ptr<RunManager> run(new RunManager(myReader, testBW, opti));

  std::cout << "Start Fit" << std::endl;
  run->startFit(par);

  std::cout << "Minimized final par :\t" << std::endl;
  std::cout << "final M:\t" << par.GetParameterValue(0) << std::endl;
  std::cout << "final T:\t" << par.GetParameterValue(1) << std::endl;

  return 0;
}
