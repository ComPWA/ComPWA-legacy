// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Test-Application of the Minuit2 Optimizer-IF.
/*!
 * @file MinuitTestApp.cpp
 * This tiny application tests the interface to the Minuit2 Optimizer. The test
 * dataset is generated in the PolyFit.hpp class, which creates smeared 1-dim data
 * according to a polynomial function. Then the Minuit2-IF is used to fit the same
 * polynomial to the smeared points and as a result the optimized parameters are
 * printed.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

// Boost header files go here
#include <boost/lexical_cast.hpp>

//#include "ErrLogger/ErrLogger.hh"

// Minimizer Interface header files go here
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"

// The toy-data to fit to
#include "executables/test/PolyFit.hpp"

using namespace std;

using ComPWA::ControlParameter;
using ComPWA::Optimizer::Optimizer;
using ComPWA::ParameterList;
using ComPWA::DoubleParameter;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt" << std::endl;
  std::cout << std::endl;

  double p0=-10., p1=10., p2=1., p3=-0.01, sigma_smear=3;

  // Generate data distribution
  //shared_ptr<ControlParameter> myFit(new PolyFit(p0, p1, p2, p3, sigma_smear));
  std::shared_ptr<ControlParameter> myFit = PolyFit::createInstance(p0, p1, p2, p3, sigma_smear);

  //--------------------------Minimizer IF --------------------------------------------------------
  vector<shared_ptr<Optimizer> > myMinimizerList;

  // Initiate parameters
  ParameterList par;
  par.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("p0",-50,-100,-5,50)));
  par.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("p1",50,0,100,50)));
  par.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("p2",10,-20,20,10)));
  par.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("p3",-0.1,-0.2,0,0.05)));

  // Add minimizers
  myMinimizerList.push_back(shared_ptr<Optimizer> (new ComPWA::Optimizer::Minuit2::MinuitIF(myFit,par)));

  // Loop over minimizers (at the moment this means: Geneva, MinuitIF or Geneva then MinuitIF)
  for(unsigned int Nmin=0; Nmin<myMinimizerList.size(); Nmin++){
    // Pointer to one ot the used minimizers
    shared_ptr<Optimizer> minimizer = myMinimizerList[Nmin];
    // Do the actual minimization
    std::shared_ptr<ComPWA::FitResult> genResult = minimizer->exec(par);

    std::cout << "Minimizer " << Nmin << "\t final par :\t" << genResult << std::endl;
    for(unsigned int i=0; i<par.GetNDouble(); i++)
    	std::cout << "final par "<< i << ":\t" << std::setprecision(9)
    << par.GetDoubleParameterValue(i) << std::endl;
    std::cout << "Done ..." << std::endl << std::endl;
  }

  // Plot results
  //myFit->drawGraph(val[0],val[1],val[2],val[3]);
  return 0;
}
