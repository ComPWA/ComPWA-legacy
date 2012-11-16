//! Test-Application of the Optimizer-IF.
/*!
 * @file OptimizerTestApp.cpp
 * This tiny application tests the interface to the Optimizers Minuit2 and Geneva.
 * The test dataset is generated in the PolyFit.hpp class, which creates smeared
 * 1-dim data according to a polynomial function. Then the Optimizer-IF implemen-
 * tations Minuit2 (OIFMinuit.hpp) and Geneva (OIFGeneva.hpp) are used one after
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
#include "OIFMinuit.hpp"
#include "OIFGeneva.hpp"
#include "PWAParameter.hpp"

// The toy-data to fit to
#include "PolyFit.hpp"

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::string whichMinimizer("all");
  double p0=-10., p1=10., p2=1., p3=-0.01, sigma_smear=3;

  // Generate data distribution
  std::shared_ptr<OIFData> myFit(new PolyFit(p0, p1, p2, p3, sigma_smear));

  //--------------------------Minimizer IF --------------------------------------------------------
  std::vector<shared_ptr<OIFBase> > myMinimizerList;

  // Add minimizers
  if (whichMinimizer=="Geneva") myMinimizerList.push_back(shared_ptr<OIFBase> (new OIFGeneva(myFit)));
  else if (whichMinimizer=="Minuit") myMinimizerList.push_back(shared_ptr<OIFBase> (new OIFMinuit(myFit)));
  else if (whichMinimizer=="all") {
    myMinimizerList.push_back(shared_ptr<OIFBase> (new OIFGeneva(myFit)));
    myMinimizerList.push_back(shared_ptr<OIFBase> (new OIFMinuit(myFit)));
  }else{
   std::cout << "Minimizer/t" << whichMinimizer << "\tdoesn't exist" << std::endl;
   return 0;
  }

  // Initiate parameters
  /*double val[4], min[4], max[4], err[4];
  val[0] = -11; max[0] = 0; min[0] = -20; err[0] = 3;
  val[1] = 9.8; max[1] = 15; min[1] = 5; err[1] = 2;
  val[2] = 1.1; max[2] = 1.5; min[2] = 0.5; err[2] = 0.3;
  val[3] = -0.008; max[3] = 0.; min[3] = -0.02; err[3] = 0.005;*/
  std::vector<PWAParameter<double> > par;
  par.push_back(PWAParameter<double>(-11,-20,0,3));
  par.push_back(PWAParameter<double>(9.8,5,15,2));
  par.push_back(PWAParameter<double>(1.1,0.5,1.5,0.3));
  par.push_back(PWAParameter<double>(-0.008,-0.02,0,0.005));

  // Loop over minimizers (at the moment this means: Geneva, Minuit or Geneva then Minuit)
  for(unsigned int Nmin=0; Nmin<myMinimizerList.size(); Nmin++){
    // Pointer to one ot the used minimizers
      std::shared_ptr<OIFBase> minimizer = myMinimizerList[Nmin];
    // Do the actual minimization
    double genResult = minimizer->exec(par);

    std::cout << "Minimizer " << Nmin << "\t final par :\t" << genResult << std::endl;
    for(unsigned int i=0; i<par.size(); i++)
      std::cout << "final par "<< i << ":\t" << par[i] << std::endl;
    std::cout << "Done ..." << std::endl << std::endl;
  }

  // Plot results
  //myFit->drawGraph(val[0],val[1],val[2],val[3]);
  return 0;
}
