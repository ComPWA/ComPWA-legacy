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
//! Test-Application of the Optimizer-IF.
/*!
 * @file OptimizerTestApp.cpp
 * This tiny application tests the interface to the Optimizers Minuit2 and Geneva.
 * The test dataset is generated in the PolyFit.hpp class, which creates smeared
 * 1-dim data according to a polynomial function. Then the Optimizer-IF implemen-
 * tations Minuit2 (MinuitIF.hpp) and Geneva (GenevaIF.hpp) are used one after
 * the other to fit the same polynomial to the smeared points. As a result the
 * optimized parameters are printed. Note: In this example Minuit2 uses the final
 * parameters of Geneva as starting values!
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
#include "Optimizer/Geneva/GenevaIF.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"

// The toy-data to fit to
#include "test/PolyFit.hpp"

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt" << std::endl;
  std::cout << std::endl;

  std::string whichMinimizer("all");
  double p0=-10., p1=10., p2=1., p3=-0.01, sigma_smear=3;
  // Generate data distribution
  std::shared_ptr<ControlParameter> myFit = PolyFit::createInstance(p0, p1, p2, p3, sigma_smear);


  ParameterList par;
  par.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("p0",-50,-100,-5,50)));
  par.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("p1",50,0,100,50)));
  par.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("p2",10,-20,20,10)));
  par.AddParameter(std::shared_ptr<DoubleParameter>(new DoubleParameter("p3",-0.1,-0.2,0,0.05)));

  //--------------------------Minimizer IF --------------------------------------------------------
  std::vector<std::shared_ptr<Optimizer> > myMinimizerList;
  // Add minimizers
  if (whichMinimizer=="Geneva") myMinimizerList.push_back(std::shared_ptr<Optimizer> (new GenevaIF(myFit)));
  else if (whichMinimizer=="Minuit") myMinimizerList.push_back(std::shared_ptr<Optimizer> (new MinuitIF(myFit,par)));
  else if (whichMinimizer=="all") {
    myMinimizerList.push_back(std::shared_ptr<Optimizer> (new GenevaIF(myFit)));
    myMinimizerList.push_back(std::shared_ptr<Optimizer> (new MinuitIF(myFit,par)));
  }else{
   std::cout << "Minimizer\t" << whichMinimizer << "\tdoesn't exist" << std::endl;
   return 0;
  }

  std::cout << "Starting Parameters:" << std::endl;
  for(unsigned int i=0; i<par.GetNDouble(); i++)
    std::cout << "Parameter "<< i << ":\t" << par.GetParameterValue(i) << std::endl;
  std::cout << std::endl << std::endl << std::endl;

  // Loop over minimizers (at the moment this means: Geneva, MinuitIF or Geneva then MinuitIF)
  for(unsigned int Nmin=0; Nmin<myMinimizerList.size(); Nmin++){
    // Pointer to one of the used minimizers
    std::shared_ptr<Optimizer> minimizer = myMinimizerList.at(Nmin);
    // Do the actual minimization
    double genResult = minimizer->exec(par);

    std::cout << "Minimizer " << Nmin << "\t final par :\t" << genResult << std::endl;

    //for(unsigned int i=0; i<par.GetNDouble(); i++)
    //  std::cout << "final par "<< i << ":\t" << par.GetParameterValue(i) << std::endl;
    std::cout << "Final: " << par << std::endl;
    std::cout << "Done ..." << std::endl << std::endl;
  }

  // Plot results
  //myFit->drawGraph(par.GetParameterValue(0),par.GetParameterValue(1),par.GetParameterValue(2),par.GetParameterValue(3));
  return 0;
}
