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
//! Test-Application for Black-Scholes stock model.
/*!
 * @file BlackScholesApp.cpp
 * This tiny application tests a simple fit procedure on stock data using
 * the Black-Scholes model.
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
#include "TGraph.h"
#include "TFile.h"

//Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"

// ComPWA header files go here
#include "DataReader/StockReader/StockReader.hpp"
#include "Physics/BlackScholes/BlackScholes.hpp"
#include "Estimator/StockLeastSquares/StockLeastSquares.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

const unsigned int nDays=800;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv){
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details see license.txt" << std::endl;
  std::cout << std::endl;

  std::string file="test/stocksHDD.csv";
  std::cout << "Load Modules" << std::endl;
  std::shared_ptr<Data> myReader(new StockReader(file ,1, nDays));
  std::shared_ptr<Amplitude> stockAmp(new BlackScholes(nDays));
 /* std::shared_ptr<ControlParameter> testEsti = StockLeastSquares::createInstance(stockAmp, myReader,0,0);
  std::shared_ptr<Optimizer> opti(new MinuitIF(testEsti));

  // Initiate parameters
  ParameterList par;
  stockAmp->fillStartParVec(par);

  std::cout << "Events :\t" << myReader->getNEvents() << std::endl << std::endl;
  std::cout << "Inital par :\t" << par << std::endl << std::endl;
  std::cout << "Inital Fitness :\t" << testEsti->controlParameter(par) << std::endl << std::endl;

  std::cout << "Start Fit" << std::endl;
  double genResult = opti->exec(par);

  std::cout << "Minimized final par :\t" << genResult << std::endl << par << std::endl << std::endl;

  //Create some output
  TGraph* sv = new TGraph(nDays);
  //sv->GetXaxis()->SetTitle("Day");
  //sv->GetXaxis()->CenterTitle();
  //sv->GetYaxis()->SetTitle("Stock Value");
  //sv->GetYaxis()->CenterTitle();

  for(unsigned int i = 0; i < myReader->getNEvents(); i++){
      Event event(myReader->getEvent(i));
      const Particle &a(event.getParticle(0));

      sv->SetPoint(i, i, a.E);
  }

  //BreitWigner *drawBW = (BreitWigner*) (&(*testBW));
  TF1* fitresult = new TF1("fitresult", ((BlackScholes*)stockAmp.get()), &BlackScholes::drawInt,0.,nDays,4,"BlackScholes","stockvalue");
  for(unsigned int i=0; i<par.GetNParameter(); i++)
    fitresult->FixParameter(i, par.GetDoubleParameter(i).GetValue());
  //fitresult->FixParameter(2, par[2]->GetValue());
  sv->Fit(fitresult);

  std::cout << "Writing output" << std::endl;

  TFile output("test/StockValueTest.root","RECREATE","ROOT_Tree");
  sv->Write();
  output.Write();
  output.Close();
*/
  return 0;
}
