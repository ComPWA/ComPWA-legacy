//! Test-Application for full fit with simple BW-dalitz-model.
/*!
 * @file DalitzFitApp.cpp
 * This tiny application tests a dalitz-fit procedure with a simple resonance
 * model. It uses the simple LH-estimator MinLogLH, it reads data
 * via the root-reader module RootReader and uses the intensity provided by
 * the Breit-Wigner-Sum  physics module AmplitudeSum. The optimization of the
 * parameters is done with the Minuit2 module MinuitIF. As result the
 * optimized parameters are printed to the terminal.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

// Root header files go here
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH2D.h"

//Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Analysis/PlotData.hpp"

const Double_t M = 3.096916; // GeV/c² (J/psi+)
const Double_t Br = 0.000093; // GeV/c² (width)
const Double_t m1 = 0.; // GeV/c² (gamma)
const Double_t m2 = 0.139570; // GeV/c² (pi)
const Double_t m3 = 0.139570; // GeV/c² (pi)
//const Double_t c = 299792458.; // m/s
const Double_t PI = 3.14159; // m/s

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::string file="gen-out.root";
  std::string outFile="fit-out.root";
  AmplitudeSetup ini("Analysis/DKsKKRes.xml");//put start parameters here
  std::cout << "Load Modules" << std::endl;
  std::shared_ptr<Data> myReader(new RootReader(file, false,"data"));
  std::shared_ptr<Data> myPHSPReader(new RootReader(file, false,"mc"));
  std::shared_ptr<Amplitude> amps(new AmpSumIntensity(M, Br, m1, m2, m3, ini));
  std::shared_ptr<ControlParameter> esti = MinLogLH::createInstance(amps, myReader, myPHSPReader);
  //std::shared_ptr<Estimator> esti(new MinLogLH(amps, myReader, myPHSPReader));
  std::shared_ptr<Optimizer> opti(new MinuitIF(esti));

  // Initiate parameters
  ParameterList par;
  amps->fillStartParVec(par); //perfect startvalues
  amps->printAmps();
  amps->getMaxVal();
  std::cout << "LH with start parameters from xml file: " << esti->controlParameter(par) << std::endl;
//  for(unsigned int i=0; i<par.GetNDouble(); i++){
//    par.GetDoubleParameter(i).SetValue(300./(double)(i+1));
//    par.GetDoubleParameter(i).SetError(10.);
//  }
//  std::cout << "LH mit folgenden intensitaeten: " << esti->controlParameter(par) << std::endl;
//
//  for(unsigned int i=0; i<par.GetNDouble(); i++){
//      std::cout << "Parameter " << i << " = " << par.GetDoubleParameter(i).GetValue() << std::endl;
//  }

 // std::cout << "Fixing 5 of 7 parameters " << std::endl;
  //for(unsigned int i=2; i<par.GetNDouble(); i++){
  //    par.GetDoubleParameter(i).FixParameter(true);
  //  }

  std::cout << "Start Fit" << std::endl;
  double genResult = opti->exec(par);
  std::cout << "Final LH = " << genResult << std::endl;
  amps->printAmps();

  plotData pl("muh",outFile, myReader,myPHSPReader,amps);
  pl.setPar(par);
  pl.plot();

  cout<<"END"<<endl;
  return 0;
}
