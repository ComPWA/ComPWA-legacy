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
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <memory>
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"
#include "Optimizer/Minuit2/MinuitIF.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"

using namespace ROOT::Minuit2;

MinuitIF::MinuitIF(std::shared_ptr<ControlParameter> theData, ParameterList& par) : _myFcn(theData, par){
  //_myFcn = new MIMinuitFcn(theData);
}

MinuitIF::~MinuitIF(){
  //std::cout << "MinuitIF::~MinuitIF: I'll be back" << std::endl;
  //delete _myFcn;
}

const double MinuitIF::exec(ParameterList& par){
  //std::string s;
  //std::stringstream out;

  MnUserParameters upar;
  for(unsigned int i=0; i<par.GetNDouble(); ++i){ //only doubles for minuit

    //out.str("");
    //out << i;
    //s = out.str();

    //use as much information as possible: (just bounds but no error not supported by minuit)
    //try{
      std::shared_ptr<DoubleParameter> actPat = par.GetDoubleParameter(i);
    //}
    if( actPat->UseBounds() && actPat->HasError() )
      upar.Add(actPat->GetName(), actPat->GetValue(), actPat->GetError(), actPat->GetMaxValue(), actPat->GetMinValue());
    else if( actPat->HasError() )
      upar.Add(actPat->GetName(), actPat->GetValue(), actPat->GetError());
    else
      upar.Add(actPat->GetName(), actPat->GetValue());

    if(actPat->IsFixed())
      upar.Fix(actPat->GetName());
  }

  MnMigrad migrad(_myFcn, upar);
  std::cout <<"start migrad "<< std::endl;
  //for(unsigned int i=0; i<par.GetNDouble(); i++)
 //   std::cout << upar.Parameter(i).Value() << " " << upar.Parameter(i).IsFixed() << std::endl;
  FunctionMinimum minMin = migrad(200,0.001);//TODO

 //if(!minMin.IsValid()) {
   //try with higher strategy
 //    std::cout <<"FM is invalid, try with strategy = 2."<< std::endl;
//   MnMigrad migrad2(_myFcn, minMin.UserState(), MnStrategy(2));
//   minMin = migrad2(10,0.1);//TODO
// }

  //save minimzed values
  for(unsigned int i=0; i<par.GetNDouble(); ++i){
    //out.str("");
    //out << i;
   // s = out.str();
    std::shared_ptr<DoubleParameter> actPat = par.GetDoubleParameter(i);
    if(!actPat->IsFixed()){
      actPat->SetValue(minMin.UserState().Value(actPat->GetName()));
      actPat->SetError(minMin.UserState().Error(actPat->GetName()));
    }
  }

  return minMin.Fval();
}

