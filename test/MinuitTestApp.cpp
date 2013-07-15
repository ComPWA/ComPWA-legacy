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
#include "PolyFit.hpp"

using namespace std;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  double p0=-10., p1=10., p2=1., p3=-0.01, sigma_smear=3;

  // Generate data distribution
  //shared_ptr<ControlParameter> myFit(new PolyFit(p0, p1, p2, p3, sigma_smear));
  shared_ptr<ControlParameter> myFit = PolyFit::createInstance(p0, p1, p2, p3, sigma_smear);

  //--------------------------Minimizer IF --------------------------------------------------------
  vector<shared_ptr<Optimizer> > myMinimizerList;

  // Add minimizers
  myMinimizerList.push_back(shared_ptr<Optimizer> (new MinuitIF(myFit)));

  // Initiate parameters
  ParameterList par;
  par.AddParameter(DoubleParameter(-11,-20,0,3));
  par.AddParameter(DoubleParameter(9.8,5,15,2));
  par.AddParameter(DoubleParameter(1.1,0.5,1.5,0.3));
  par.AddParameter(DoubleParameter(-0.008,-0.02,0,0.005));

  // Loop over minimizers (at the moment this means: Geneva, MinuitIF or Geneva then MinuitIF)
  for(unsigned int Nmin=0; Nmin<myMinimizerList.size(); Nmin++){
    // Pointer to one ot the used minimizers
    shared_ptr<Optimizer> minimizer = myMinimizerList[Nmin];
    // Do the actual minimization
    double genResult = minimizer->exec(par);

    std::cout << "Minimizer " << Nmin << "\t final par :\t" << genResult << std::endl;
    for(unsigned int i=0; i<par.GetNDouble(); i++)
      std::cout << "final par "<< i << ":\t" << par.GetParameterValue(i) << std::endl;
    std::cout << "Done ..." << std::endl << std::endl;
  }

  // Plot results
  //myFit->drawGraph(val[0],val[1],val[2],val[3]);
  return 0;
}
