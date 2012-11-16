//! Test-Application for full fit with simple modules.
/*!
 * @file IntegrationTestApp.cpp
 * This tiny application tests a simple fit procedure with a set of simple
 * modules. It uses a simle \f$\chi^{2}\f$ estimator EIFChiOneD, it reads data
 * via the root-reader module DIFRootReader and uses the intensity provided by
 * the simple 1D-Breit-Wigner physics module PIFBW. The optimization of the
 * parameters is done with the Minuit2 module OIFMinuit. As result the
 * optimized parameters are printed to the terminal.
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
#include "OIFMinuit.hpp"
#include "PWAParameter.hpp"

//Test header files go here
#include "PolyFit.hpp"

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  string file="test/2Part-4vecs.root";
  std::cout << "Load Modules" << std::endl;
  std::shared_ptr<DIFBase> myReader(new DIFRootReader(file));
  std::shared_ptr<PIFBase> testBW(new PIFBW());
  std::shared_ptr<EIFBase> testEsti(new EIFChiOneD(testBW, myReader));
  std::shared_ptr<OIFBase> opti(new OIFMinuit(testEsti));

  // Initiate parameters
  //double val[2], min[2], max[2], err[2];
  //val[0] = 1.5; max[0] = 2.5; min[0] = 0.5; err[0] = 0.5;
  //val[1] = 0.3; max[1] = 0.5; min[1] = 0.1; err[1] = 0.1;
  std::vector<PWAParameter<double> > par;
  par.push_back(PWAParameter<double>(1.5,0.5,2.5,0.5));
  par.push_back(PWAParameter<double>(0.3,0.1,0.5,0.1));

  std::cout << "Start Fit" << std::endl;
  double genResult = opti->exec(par);

  std::cout << "Minimized final par :\t" << genResult << std::endl;
  std::cout << "final M:\t" << par[0].GetValue() << " +- " << par[0].GetError() << std::endl;
  std::cout << "final T:\t" << par[1].GetValue() << " +- " << par[1].GetError() << std::endl;

  return 0;
}