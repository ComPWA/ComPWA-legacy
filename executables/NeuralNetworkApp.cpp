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
#include "Physics/NeuralNetwork/NeuralNetwork.hpp"
#include "Estimator/StockLeastSquares/StockLeastSquares.hpp"
#include "Optimizer/Geneva/GenevaIF.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

//const unsigned int nDays=13;
const unsigned int nInputDays=50;
const unsigned int nTestDays=5;
const unsigned int nPredDays=2;

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
  std::shared_ptr<Data> myReader(new StockReader(file ,1, nInputDays+nTestDays));
  std::shared_ptr<NeuralNetwork> stockAmp(new NeuralNetwork(nInputDays, nTestDays,3,0.025));
  // Initiate parameters
  ParameterList par;
  stockAmp->fillStartParVec(par);
  std::shared_ptr<ControlParameter> testEsti = StockLeastSquares::createInstance(stockAmp, myReader);
  std::shared_ptr<Optimizer> opti(new GenevaIF(testEsti));
  std::shared_ptr<Optimizer> optiAlt(new MinuitIF(testEsti,par));

  std::cout << "Events :\t" << myReader->getNEvents() << std::endl << std::endl;
  if(nInputDays<11)
    std::cout << "Inital par :\t" << par << std::endl << std::endl;
  else
    std::cout << "Num Inital par :\t" << par.GetNDouble() << std::endl << std::endl;
  std::vector<double> x;
  for(unsigned int i=0; i<nInputDays; i++) x.push_back(i);
  std::cout << "Inital OutPar (input =1):\t" << stockAmp->intensity(x,par).GetParameterValue(0) << std::endl << std::endl;
  std::cout << "Inital Fitness :\t" << testEsti->controlParameter(par) << std::endl << std::endl;

  std::cout << "Start Fit" << std::endl;
  double genResult = opti->exec(par);

  std::cout << "Minimized final par :\t" << genResult << std::endl << std::endl;

  //-----------------------------Create some output
  TGraph* sv = new TGraph(nTestDays);
  sv->SetName("TestValues");
  sv->SetTitle("TestValues");
  TGraph* sv2 = new TGraph(nTestDays);
  sv2->SetName("NNValues");
  sv2->SetTitle("NNValues");
  TGraph* delta = new TGraph(nTestDays);
  delta->SetName("NN_Test_Diff");
  delta->SetTitle("NN_Test_Diff");

  std::vector<double> d;
  for(unsigned int k=0; k<nInputDays; k++){
      Event event(myReader->getEvent(k));
      const Particle &a(event.getParticle(0));
      d.push_back(a.E);
  }
  ParameterList outVals = stockAmp->intensity(d,par);
  for(unsigned int j = 0; j < nTestDays; j++){
      Event event(myReader->getEvent(j+nInputDays));
      const Particle &a(event.getParticle(0));

      std::cout << "Day " << j+nInputDays << "\t SV " << a.E << "\t NNV " << outVals.GetParameterValue(j) << " " << outVals.GetDoubleParameter(j)->GetName() << std::endl;

      sv->SetPoint(j, j+nInputDays, a.E);

      sv2->SetPoint(j, j+nInputDays, outVals.GetParameterValue(j));

      delta->SetPoint(j, j+nInputDays, (a.E-outVals.GetParameterValue(j)));
  }

  //------------------------------Predict some values
  TGraph* futureVal = new TGraph(nPredDays);
  futureVal->SetName("Future Values");
  futureVal->SetTitle("Future Values");
  TGraph* predVal = new TGraph(nPredDays);
  predVal->SetName("Predicted Values");
  predVal->SetTitle("Predicted Values");
  TGraph* deltaPred = new TGraph(nPredDays);
  deltaPred->SetName("Futute_Pred_Diff");
  deltaPred->SetTitle("Future_Pred_Diff");

  std::shared_ptr<Data> myReader2(new StockReader(file ,1, nInputDays+nTestDays+nPredDays));
  std::vector<double> f;
  for(unsigned int k=nPredDays; k<(nInputDays+nPredDays); k++){
      Event event(myReader2->getEvent(k));
      const Particle &a(event.getParticle(0));
      f.push_back(a.E);
  }
  ParameterList predVals = stockAmp->intensity(f,par);
  for(unsigned int j = 0; j < nPredDays; j++){
      Event event(myReader2->getEvent(j+nInputDays+nTestDays));
      const Particle &a(event.getParticle(0));

      std::cout << "Future Day " << j+nInputDays+nTestDays << "\t SV " << a.E << "\t NNV " << predVals.GetParameterValue(j) << " " << predVals.GetDoubleParameter(j)->GetName() << std::endl;

      futureVal->SetPoint(j, j+nInputDays+nTestDays, a.E);

      predVal->SetPoint(j, j+nInputDays+nTestDays, predVals.GetParameterValue(j));

      deltaPred->SetPoint(j, j+nInputDays+nTestDays, (a.E-predVals.GetParameterValue(j)));
  }


  std::cout << "Writing output" << std::endl;

  TFile output("test/NeuralNetworkTest.root","RECREATE","ROOT_Tree");
  sv->Write();
  sv2->Write();
  delta->Write();
  futureVal->Write();
  predVal->Write();
  deltaPred->Write();
  output.Write();
  output.Close();

  return 0;
}
