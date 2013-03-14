//! Test-Application for full fit with simple modules.
/*!
 * @file IntegrationTestApp.cpp
 * This tiny application tests a simple fit procedure with a set of simple
 * modules. It uses a simle \f$\chi^{2}\f$ estimator ChiOneD, it reads data
 * via the root-reader module RootReader and uses the intensity provided by
 * the simple 1D-Breit-Wigner physics module BreitWigner. The optimization of the
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

//Core header files go here
#include "Core/PWAEvent.hpp"
#include "Core/PWAParticle.hpp"
#include "Core/PWAParameter.hpp"
#include "Core/PWAGenericPar.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/BreitWigner/BreitWigner.hpp"
#include "Estimator/ChiOneD/ChiOneD.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::string file="test/2Part-4vecs.root";
  std::cout << "Load Modules" << std::endl;
  std::shared_ptr<Data> myReader(new RootReader(file, true));
  std::shared_ptr<Amplitude> testBW(new BreitWigner(0.,5.));
  std::shared_ptr<Estimator> testEsti(new ChiOneD(testBW, myReader));
  std::shared_ptr<Optimizer> opti(new MinuitIF(testEsti));

  // Initiate parameters
  ParameterList par;
  testBW->fillStartParVec(par);
  par.GetDoubleParameter(0).SetValue(1.7);
  par.GetDoubleParameter(0).SetValue(0.2);

  std::cout << "Inital par :\t" << std::endl;
  std::cout << "inital M:\t" << par.GetDoubleParameter(0).GetValue() << std::endl;
  std::cout << "inital T:\t" << par.GetDoubleParameter(1).GetValue() << std::endl;

  std::cout << "Start Fit" << std::endl;
  double genResult = opti->exec(par);

  std::cout << "Minimized final par :\t" << genResult << std::endl;
  std::cout << "final M:\t" << par.GetDoubleParameter(0).GetValue() << std::endl;
  std::cout << "final T:\t" << par.GetDoubleParameter(1).GetValue() << std::endl;

  //Create some output
  TH1D* bw = new TH1D("bw","inv. mass of 2 particles",1000,0.,2.4);
  bw->GetXaxis()->SetTitle("m_{12} / GeV");
  bw->GetXaxis()->CenterTitle();
  bw->GetYaxis()->SetTitle("#");
  bw->GetYaxis()->CenterTitle();

  for(unsigned int i = 0; i < myReader->getNEvents(); i++){
      PWAEvent event;
      PWAParticle a, b;
      double masssq = 0;

      if(!myReader->getEvent(i, event)) continue;
      //if(!event) continue;
      //cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles() << endl;
      event.getParticle(0,a);
      event.getParticle(1,b);
      masssq = pow(a.getE()+b.getE(),2) - pow(a.getPx()+b.getPx() ,2) - pow(a.getPy()+b.getPy() ,2) - pow(a.getPz()+b.getPz() ,2);

      bw->Fill(sqrt(masssq));
  }

  //BreitWigner *drawBW = (BreitWigner*) (&(*testBW));
  TF1* fitresult = new TF1("fitresult", ((BreitWigner*)testBW.get()), &BreitWigner::drawInt,0.,2.4,3,"PIFBW","intensity");
  fitresult->FixParameter(0, par.GetDoubleParameter(0).GetValue());
  fitresult->FixParameter(1, par.GetDoubleParameter(1).GetValue());
  //fitresult->FixParameter(2, par[2]->GetValue());
  bw->Fit(fitresult);

  TFile output("test/IntegrationTestBinned.root","RECREATE","ROOT_Tree");
  bw->Write();
  output.Write();
  output.Close();

  return 0;
}
