//! Test-Application for fit with simple modules using RunManager.
/*!
 * @file RunManagerTestApp.cpp
 * This tiny application tests a simple fit procedure with a set of simple
 * modules. It uses a simle \f$\chi^{2}\f$ estimator EIFChiOneD, it reads data
 * via the root-reader module DIFRootReader and uses the intensity provided by
 * the simple 1D-Breit-Wigner physics module PIFBW. The optimization of the
 * parameters is done with the Minuit2 module OIFMinuit. The class RunManager
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
#include "DIFRootReader.hpp"
#include "PIFBW.hpp"
#include "EIFChiOneD.hpp"
#include "OIFMinuit.hpp"
#include "PWAParameter.hpp"
#include "RunManager.hpp"

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
  std::shared_ptr<EIFBase> testEsti(new EIFChiOneD(testBW, myReader)); //TODO: <- should be done by runManager
  std::shared_ptr<OIFBase> opti(new OIFMinuit(testEsti));
  std::shared_ptr<RunManager> run(new RunManager(myReader, testEsti, testBW, opti));

  // Initiate parameters
  std::vector<PWAParameter<double> > par;
  par.push_back(PWAParameter<double>(1.5,0.5,2.5,0.5));
  par.push_back(PWAParameter<double>(0.3,0.1,0.5,0.1));

  std::cout << "Start Fit" << std::endl;
  run->startFit(par);

  std::cout << "Minimized final par :\t" << std::endl;
  std::cout << "final M:\t" << par[0].GetValue() << " +- " << par[0].GetError() << std::endl;
  std::cout << "final T:\t" << par[1].GetValue() << " +- " << par[1].GetError() << std::endl;

  return 0;
}
