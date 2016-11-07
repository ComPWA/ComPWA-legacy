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
//! Test-Application for full fit with simple modules.
/*!
 * @file IntegrationTestApp.cpp
 * This tiny application tests a simple fit procedure with a set of simple
 * modules. It uses a simle LH-estimator MinLogLH, it reads data
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
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

// ComPWA header files go here
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/BreitWigner/BreitWigner.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

using namespace ComPWA;
using namespace DataReader;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt" << std::endl;
  std::cout << std::endl;

  std::string file="test/2Part-4vecs.root";
  std::cout << "Load Modules" << std::endl;
  std::shared_ptr<Data> myReader(new RootReader(file, "data"));
  std::shared_ptr<Amplitude> testBW(
		  new Physics::BreitWigner::BreitWigner(0.,5.)
  );

  // Initiate parameters
  ParameterList par;
  testBW->FillParameterList(par);

  std::shared_ptr<Optimizer::ControlParameter> testEsti =
		  Estimator::MinLogLH::MinLogLH::createInstance(
				  testBW, myReader, std::shared_ptr<Data>(),
				  std::shared_ptr<Data>()
		  );
  std::shared_ptr<Optimizer::Optimizer> opti(
		  new Optimizer::Minuit2::MinuitIF(testEsti, par)
  );

  par.GetDoubleParameter(0)->SetValue(1.7);
  par.GetDoubleParameter(0)->SetValue(0.2);

  std::cout << "Inital par :\t" << std::endl;
  std::cout << "inital M:\t" << par.GetDoubleParameter(0)->GetValue() << std::endl;
  std::cout << "inital T:\t" << par.GetDoubleParameter(1)->GetValue() << std::endl;
  std::cout << "Start Fit" << std::endl;
  std::shared_ptr<FitResult> genResult = opti->exec(par);

  std::cout << "Minimized final par :\t" << genResult->getResult()<< std::endl;
  std::cout << "final M:\t" << par.GetDoubleParameter(0)->GetValue() << std::endl;
  std::cout << "final T:\t" << par.GetDoubleParameter(1)->GetValue() << std::endl;

  //Create some output
  TH1D* bw = new TH1D("bw","inv. mass of 2 particles",1000,0.,2.4);
  bw->GetXaxis()->SetTitle("m_{12} / GeV");
  bw->GetXaxis()->CenterTitle();
  bw->GetYaxis()->SetTitle("#");
  bw->GetYaxis()->CenterTitle();

  for(unsigned int i = 0; i < myReader->getNEvents(); i++){
      Event event(myReader->getEvent(i));
      double masssq = 0;

      //if(!myReader->getEvent(i, event)) continue; TODO: try wxception
      //if(!event) continue;
      //cout << "Event: \t" << i << "\t NParticles: \t" << event.getNParticles() << endl;
      const Particle &a(event.getParticle(0));
      const Particle &b(event.getParticle(1));
      masssq = pow(a.E+b.E,2) - pow(a.px+b.px ,2) - pow(a.py+b.py ,2) - pow(a.pz+b.pz ,2);

      bw->Fill(sqrt(masssq));
  }

  //BreitWigner *drawBW = (BreitWigner*) (&(*testBW));
  TF1* fitresult = new TF1(
		  "fitresult",
		  ((BreitWigner*) testBW.get()),
		  &Physics::BreitWigner::BreitWigner::drawInt,
		  0.,2.4,3,
		  "PIFBW",
		  "intensity"
  );

  fitresult->FixParameter(0, par.GetDoubleParameter(0)->GetValue());
  fitresult->FixParameter(1, par.GetDoubleParameter(1)->GetValue());
  //fitresult->FixParameter(2, par[2]->GetValue());
  bw->Fit(fitresult);

  TFile output("test/IntegrationTest.root","RECREATE","ROOT_Tree");
  bw->Write();
  output.Write();
  output.Close();

  return 0;
}
